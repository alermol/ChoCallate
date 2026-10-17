#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: write_versions.sh [OPTIONS] OUTPUT TOOL [TOOL ...]

Writes the versions of the requested tools to OUTPUT.

  --format yaml|flat      tool: "version" (default) or tool=version lines
  --pipeline-version VER  also record choCallate: "VER"
EOF
}

format="yaml"
pipeline_version=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --format)
            format="$2"
            shift 2
            ;;
        --pipeline-version)
            pipeline_version="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

if [[ $# -lt 2 ]]; then
    usage >&2
    exit 2
fi
if [[ "${format}" != "yaml" && "${format}" != "flat" ]]; then
    echo "--format must be yaml or flat" >&2
    exit 2
fi

output="$1"
shift

version_of() {
    local tool="$1"
    local version=""

    if ! command -v "${tool}" >/dev/null 2>&1 &&
        [[ "${tool}" != "pysam" ]]; then
        echo "write_versions.sh: '${tool}' is not on PATH" >&2
        exit 3
    fi

    case "${tool}" in
        bcftools)
            version=$(bcftools --version 2>/dev/null |
                awk 'NR==1 {print $2}')
            ;;
        samtools)
            version=$(samtools --version 2>/dev/null |
                awk 'NR==1 {print $2}')
            ;;
        bgzip)
            version=$(bgzip --version 2>/dev/null |
                awk 'NR==1 {print $3}')
            ;;
        tabix)
            version=$(tabix --version 2>/dev/null |
                awk 'NR==1 {print $3}')
            ;;
        bwa)
            version=$( { bwa 2>&1 || true; } |
                awk '/^Version:/ {print $2; exit}')
            ;;
        bowtie2)
            version=$(bowtie2 --version 2>/dev/null |
                awk '/ version / {print $NF; exit}')
            ;;
        minimap2)
            version=$(minimap2 --version 2>/dev/null |
                awk 'NR==1 {print $1}')
            ;;
        freebayes)
            version=$(freebayes --version 2>/dev/null |
                sed -n \
                    's/^version:[[:space:]]*v\{0,1\}\([0-9A-Za-z.-]*\).*/\1/p' |
                head -1)
            ;;
        gatk)
            version=$(gatk --version 2>&1 |
                sed -n 's/.*GATK) v\([0-9.]*\).*/\1/p' |
                head -1)
            ;;
        mosdepth)
            version=$(mosdepth --version 2>/dev/null |
                awk 'NR==1 {print $2}')
            ;;
        bedops)
            version=$(bedops --version 2>/dev/null |
                sed -n 's/^[[:space:]]*version:[[:space:]]*\([0-9.]*\).*/\1/p' |
                head -1)
            ;;
        python)
            version=$(python --version 2>&1 |
                awk 'NR==1 {print $2}')
            ;;
        pysam)
            version=$(python -c 'import pysam; print(pysam.__version__)' \
                2>/dev/null)
            ;;
        parallel)
            version=$(parallel --version 2>/dev/null |
                awk 'NR==1 {print $3}')
            ;;
        *)
            echo "write_versions.sh: unsupported tool '${tool}'" >&2
            exit 2
            ;;
    esac

    if [[ -z "${version}" ]]; then
        echo "write_versions.sh: cannot determine version of '${tool}'" >&2
        exit 4
    fi

    if [[ "${format}" == "flat" ]]; then
        printf '%s=%s\n' "${tool}" "${version}"
    else
        printf '%s: "%s"\n' "${tool}" "${version}"
    fi
}

tmpfile=$(mktemp)
trap 'rm -f "${tmpfile}"' EXIT
for tool in "$@"; do
    version_of "${tool}" >> "${tmpfile}"
done
if [[ -n "${pipeline_version}" ]]; then
    if [[ "${format}" == "flat" ]]; then
        printf 'ChoCallate=%s\n' "${pipeline_version}" >> "${tmpfile}"
    else
        printf 'ChoCallate: "%s"\n' "${pipeline_version}" >> "${tmpfile}"
    fi
fi
sort -u "${tmpfile}" > "${output}"
