#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: postprocess_bam.sh --input BAM|-

  --output PATH          Final BAM (default: mapping.bam)
  --reference FASTA      Reference FASTA (required with --left-align)
  --dict DICT            Sequence dictionary (required with --left-align)
  --fai PATH             Reference .fai index. Verify the output BAM
                         @SQ entries against it; with --dict, compare
                         M5 checksums too.
  --mapq INT             Minimum mapping quality (default: 0)
  --threads INT          Worker threads, split across pipeline stages
                         (default: 1)
  --memory-mb INT        Task memory in MiB (sort uses half of the
                         remainder after --reserve-mb)
  --reserve-mb INT       Memory in MiB reserved for the aligner index
                         and other co-resident processes; excluded
                         before splitting sort memory (default: 0)
  --java-options OPTS    JVM options for GATK LeftAlignIndels
  --rmdup-mode MODE      Duplicate handling: none, mark or remove
                         (default: none)
  --rmdup                Alias for --rmdup-mode remove
  --left-align           Left-align indels with GATK LeftAlignIndels
EOF
}

input=""
output="mapping.bam"
reference=""
dict=""
fai=""
mapq=0
threads=1
memory_mb=0
reserve_mb=0
java_options=""
rmdup_mode="none"
left_align=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            input="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --reference)
            reference="$2"
            shift 2
            ;;
        --dict)
            dict="$2"
            shift 2
            ;;
        --fai)
            fai="$2"
            shift 2
            ;;
        --mapq)
            mapq="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --memory-mb)
            memory_mb="$2"
            shift 2
            ;;
        --reserve-mb)
            reserve_mb="$2"
            shift 2
            ;;
        --java-options)
            java_options="$2"
            shift 2
            ;;
        --rmdup-mode)
            rmdup_mode="$2"
            shift 2
            ;;
        --rmdup)
            rmdup_mode="remove"
            shift
            ;;
        --left-align)
            left_align=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$input" ]]; then
    echo "--input is required" >&2
    exit 2
fi
if [[ "$input" != "-" && ! -f "$input" ]]; then
    echo "Input BAM not found: $input" >&2
    exit 2
fi
if [[ "$left_align" == true && ( -z "$reference" || -z "$dict" ) ]]; then
    echo "--reference and --dict are required with --left-align" >&2
    exit 2
fi
if [[ -n "$fai" && ! -f "$fai" ]]; then
    echo "Reference index not found: $fai" >&2
    exit 2
fi
case "$rmdup_mode" in
    none|mark|remove) ;;
    *)
        echo "--rmdup-mode must be one of none, mark or remove" >&2
        exit 2
        ;;
esac
if [[ "$threads" -lt 1 ]]; then
    echo "--threads must be >= 1" >&2
    exit 2
fi
if [[ "$memory_mb" -lt 0 ]]; then
    echo "--memory-mb must be >= 0" >&2
    exit 2
fi
if [[ "$reserve_mb" -lt 0 ]]; then
    echo "--reserve-mb must be >= 0" >&2
    exit 2
fi
if [[ "$memory_mb" -gt 0 && "$reserve_mb" -ge "$memory_mb" ]]; then
    echo "--reserve-mb ($reserve_mb) must be smaller than" \
        "--memory-mb ($memory_mb)" >&2
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmpdir="$(mktemp -d -p "$PWD" .postprocess.XXXXXX)"
cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT

# The stages of each pipeline below run concurrently, so divide the requested
# threads between them instead of giving every stage the full count.  The
# streaming stages get one thread (or none when the budget is smaller than the
# stage count) and sort, the CPU-heavy stage, receives the remainder.
if [[ "$rmdup_mode" != none ]]; then
    streaming_stages=4    # collate, fixmate, markdup and view
else
    streaming_stages=1    # view
fi
if [[ "$threads" -gt "$streaming_stages" ]]; then
    stage_threads=1
    sort_threads=$(( threads - streaming_stages ))
else
    stage_threads=0
    sort_threads="$threads"
fi

# The aligner index is resident while sort runs in the same mapping task.
# Subtract the reserved memory first, then share half of the remainder
# between the sorting threads and leave the rest for the aligner working set
# and the other BAM-processing stages.
sort_memory_mb=768
if [[ "$memory_mb" -gt 0 ]]; then
    available_mb=$(( memory_mb - reserve_mb ))
    sort_memory_mb=$(( available_mb / 2 / sort_threads ))
    if [[ "$sort_memory_mb" -lt 1 ]]; then
        sort_memory_mb=1
    fi
fi

filtered_bam="$tmpdir/filtered.bam"
if [[ "$rmdup_mode" != none ]]; then
    markdup_flags=()
    if [[ "$rmdup_mode" == remove ]]; then
        markdup_flags=(-r)
    fi
    # Stream markdup into the MAPQ filter instead of writing dedup.bam.
    samtools collate --threads "$stage_threads" -O -u \
            -T "$tmpdir/collate" "$input" \
        | samtools fixmate --threads "$stage_threads" -m -u - - \
        | samtools sort --threads "$sort_threads" \
            -m "${sort_memory_mb}M" -u -T "$tmpdir/sort" - \
        | samtools markdup --threads "$stage_threads" -u \
            "${markdup_flags[@]}" -T "$tmpdir/markdup" - - \
        | samtools view --threads "$stage_threads" -F 4 -q "$mapq" \
            -b -o "$filtered_bam" -
else
    # Filtering before sorting reduces the number of records sort has to
    # order; no duplicate information is needed in this path.
    samtools view --threads "$stage_threads" -F 4 -q "$mapq" -u \
            "$input" \
        | samtools sort --threads "$sort_threads" \
            -m "${sort_memory_mb}M" -T "$tmpdir/sort" \
            -o "$filtered_bam" -
fi

# Fail before left-alignment and indexing when the aligner used an index
# built for another assembly.
if [[ -n "$fai" ]]; then
    check_args=("$filtered_bam" "$fai")
    if [[ -n "$dict" ]]; then
        check_args+=(--dict "$dict")
    fi
    "$script_dir/check_bam_reference.py" "${check_args[@]}"
fi

if [[ "$left_align" == true ]]; then
    rm -f "$output"
    gatk_opts=()
    if [[ -n "$java_options" ]]; then
        gatk_opts=(--java-options "$java_options")
    fi
    gatk LeftAlignIndels "${gatk_opts[@]}" \
        -I "$filtered_bam" \
        -O "$output" \
        -R "$reference" \
        --sequence-dictionary "$dict" \
        -OBI false
else
    mv "$filtered_bam" "$output"
fi

samtools index --threads "$threads" --csi "$output"
