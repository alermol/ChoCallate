#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: merge_outputs.sh --file-list FILE --reference FASTA --output PATH

  --format bcf|vcf         Output format (default: bcf)
  --cpus INT               Threads for bcftools and the chunked merge (default: 1)
  --remove-invariant       Drop monomorphic records
  --split-multiallelic     Keep one record per ALT
  --samples-dropped LIST   Comma-separated ids added to ##ChoCallate
  --bcftools PATH          bcftools executable (default: bcftools)
  --python PATH            Python interpreter for helper scripts (default: python3)
EOF
}

file_list=""
reference=""
output=""
format="bcf"
cpus=1
remove_invariant=false
split_multiallelic=false
samples_dropped=""
bcftools="${BCFTOOLS:-bcftools}"
python="${PYTHON:-python3}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --file-list)
            file_list="$2"
            shift 2
            ;;
        --reference)
            reference="$2"
            shift 2
            ;;
        --output)
            output="$2"
            shift 2
            ;;
        --format)
            format="$2"
            shift 2
            ;;
        --cpus)
            cpus="$2"
            shift 2
            ;;
        --remove-invariant)
            remove_invariant=true
            shift
            ;;
        --split-multiallelic)
            split_multiallelic=true
            shift
            ;;
        --samples-dropped)
            samples_dropped="$2"
            shift 2
            ;;
        --bcftools)
            bcftools="$2"
            shift 2
            ;;
        --python)
            python="$2"
            shift 2
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

if [[ -z "$file_list" || -z "$reference" || -z "$output" ]]; then
    echo "--file-list, --reference and --output are required" >&2
    exit 2
fi
if [[ "$format" != "bcf" && "$format" != "vcf" ]]; then
    echo "--format must be bcf or vcf" >&2
    exit 2
fi
if [[ "$cpus" -lt 1 ]]; then
    echo "--cpus must be >= 1" >&2
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmpdir="$(mktemp -d -p "${TMPDIR:-/tmp}" merge_outputs.XXXXXX)"
cleanup() {
    rm -rf "$tmpdir"
}
trap cleanup EXIT

mapfile -t files < "$file_list"
if [[ ${#files[@]} -eq 0 ]]; then
    echo "No input BCF files provided" >&2
    exit 1
fi

first="${files[0]}"

# Patch the first sample's structured header once; the merged header is
# rebuilt from it, and a single-file merge keeps it unchanged otherwise.
if [[ -n "$samples_dropped" ]]; then
    "$bcftools" view -h "$first" > "$tmpdir/first.hdr"
    awk -v dropped="$samples_dropped" \
        '/^##ChoCallate=</ {
            sub(/>$/, ",SamplesDropped=\"" dropped "\">")
        } 1' \
        "$tmpdir/first.hdr" > "$tmpdir/first.mod.hdr"
    "$bcftools" reheader -h "$tmpdir/first.mod.hdr" \
        -o "$tmpdir/first.bcf" "$first"
    files[0]="$tmpdir/first.bcf"
    printf '%s\n' "${files[@]}" > "$tmpdir/files.txt"
    file_list="$tmpdir/files.txt"
    first="${files[0]}"
fi

fill_tags="AN,AC,AF,NS,AC_Hom,AC_Het,MAF,TYPE,F_MISSING,DP:1=int(sum(FORMAT/DP))"
output_flags=(-Ob -o "$output")
if [[ "$format" == "vcf" ]]; then
    output_flags=(-Oz -o "$output")
fi

invariant_stream() {
    if [[ "$remove_invariant" == true ]]; then
        "$bcftools" filter --threads "$cpus" \
            -e 'ALT="." || ALT="<*>" || COUNT(GT="alt")=0' -Ou
    else
        "$bcftools" view --threads "$cpus" -Ou
    fi
}

join_stream() {
    if [[ "$split_multiallelic" == true ]]; then
        cat
    else
        "$bcftools" norm --threads "$cpus" -m +any -Ou
    fi
}

tag_stream() {
    "$bcftools" +fill-tags -Ou -- -t "$fill_tags"
}

annotate_stream() {
    if [[ -s "$tmpdir/cho_meta.hdr" ]]; then
        "$bcftools" annotate -h "$tmpdir/cho_meta.hdr" -Ou
    else
        cat
    fi
}

output_stream() {
    "$bcftools" view "${output_flags[@]}"
}

if [[ ${#files[@]} -eq 1 ]]; then
    "$bcftools" view --threads "$cpus" -Ou "$first" \
        | invariant_stream \
        | join_stream \
        | tag_stream \
        | output_stream
else
    for bcf in "${files[@]}"; do
        "$bcftools" index -f --csi "$bcf"
    done

    # Only the structured metadata line belongs to the merge. The command
    # line is attached to the final output by MERGE_OUTPUTS.
    "$bcftools" view -h "$first" \
        | grep '^##ChoCallate=' > "$tmpdir/cho_meta.hdr" || true

    "$python" "$script_dir/merge_bcf_tree.py" \
        --gvcf "$reference" \
        --cpus "$cpus" \
        --file-list "$file_list" \
        | "$bcftools" norm -m -any --threads "$cpus" -Ou \
        | "$python" "$script_dir/drop_refonly_duplicates.py" \
        | invariant_stream \
        | join_stream \
        | tag_stream \
        | annotate_stream \
        | output_stream
fi
