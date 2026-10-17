#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: postprocess_bam.sh --input BAM|-

  --output PATH          Final BAM (default: mapping.bam)
  --reference FASTA      Reference FASTA (required with --left-align)
  --dict DICT            Sequence dictionary (required with --left-align)
  --mapq INT             Minimum mapping quality (default: 0)
  --threads INT          Worker threads (default: 1)
  --java-options OPTS    JVM options for GATK LeftAlignIndels
  --rmdup                Remove duplicates with samtools markdup
  --left-align           Left-align indels with GATK LeftAlignIndels
EOF
}

input=""
output="mapping.bam"
reference=""
dict=""
mapq=0
threads=1
java_options=""
rmdup=false
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
        --mapq)
            mapq="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --java-options)
            java_options="$2"
            shift 2
            ;;
        --rmdup)
            rmdup=true
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
if [[ "$threads" -lt 1 ]]; then
    echo "--threads must be >= 1" >&2
    exit 2
fi

tmpdir="$(mktemp -d -p "$PWD" .postprocess.XXXXXX)"
cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT

if [[ "$rmdup" == true ]]; then
    samtools collate --threads "$threads" -O -u -T "$tmpdir/collate" "$input" \
        | samtools fixmate --threads "$threads" -m -u - - \
        | samtools sort --threads "$threads" -u -T "$tmpdir/sort" - \
        | samtools markdup --threads "$threads" -r -T "$tmpdir/markdup" - "$tmpdir/dedup.bam"
    filtered_source="$tmpdir/dedup.bam"
else
    samtools sort --threads "$threads" -o "$tmpdir/sorted.bam" "$input"
    filtered_source="$tmpdir/sorted.bam"
fi

samtools view --threads "$threads" -F 4 -q "$mapq" -b \
    -o "$tmpdir/filtered.bam" "$filtered_source"

if [[ "$left_align" == true ]]; then
    rm -f "$output"
    gatk_opts=()
    if [[ -n "$java_options" ]]; then
        gatk_opts=(--java-options "$java_options")
    fi
    gatk LeftAlignIndels "${gatk_opts[@]}" \
        -I "$tmpdir/filtered.bam" \
        -O "$output" \
        -R "$reference" \
        --sequence-dictionary "$dict" \
        -OBI false
else
    mv "$tmpdir/filtered.bam" "$output"
fi

samtools index --threads "$threads" --csi "$output"
