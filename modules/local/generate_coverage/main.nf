process GENERATE_COVERAGE {
    label 'process_low'
    cpus 1
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam"), path("tmp/input.bam.csi")
    path(include_bed)
    path(exclude_bed)

    output:
    tuple val(sample_id), path("coverage.bed"), emit: coverage

    script:
    def min_coverage = params.coverage.min_coverage != null
        ? params.coverage.min_coverage as int
        : 5
    def quantize = min_coverage > 0 ? "0:${min_coverage}:" : "0:1:"
    def keep_bin = min_coverage > 0 ? "${min_coverage}:inf" : "1:inf"
    def on_empty = params.coverage.on_empty ?: 'fail'
    def empty_message = "Sample ${sample_id}: coverage.bed is empty: " +
        "no region reaches coverage.min_coverage=${min_coverage}"
    def sort_include = include_bed.name != "N1_FILE" ? "sort-bed \"${include_bed}\" > tmp/include.sorted.bed" : ""
    def sort_exclude = exclude_bed.name != "NO_FILE" ? "sort-bed \"${exclude_bed}\" > tmp/exclude.sorted.bed" : ""
    """
    set -o pipefail
    mkdir -p tmp

    ${sort_include}
    ${sort_exclude}

    # mosdepth emits interval runs with a quantized depth label.  Keep the run
    # whose lower bound is coverage.min_coverage, then restrict to the requested
    # regions.  --no-per-base avoids the per-base intermediate file.
    mosdepth \
        --no-per-base \
        --quantize ${quantize} \
        --threads ${task.cpus} \
        tmp/mosdepth \
        tmp/input.bam

    zcat tmp/mosdepth.quantized.bed.gz \\
        | awk -v OFS='\\t' -v keep="${keep_bin}" '\$4 == keep {print \$1, \$2, \$3}' \\
        > tmp/quantized.bed
    sort-bed tmp/quantized.bed > tmp/quantized.sorted.bed

    if [ -f tmp/include.sorted.bed ]; then
        bedops --intersect tmp/quantized.sorted.bed tmp/include.sorted.bed > tmp/selected.bed
    else
        cp tmp/quantized.sorted.bed tmp/selected.bed
    fi

    if [ -f tmp/exclude.sorted.bed ]; then
        bedops --difference tmp/selected.bed tmp/exclude.sorted.bed > tmp/filtered.bed
    else
        cp tmp/selected.bed tmp/filtered.bed
    fi

    bedops --merge tmp/filtered.bed > coverage.bed

    # bedops writes an empty file when no run reaches the threshold.  Fail with
    # the sample id and threshold, or let main.nf drop the sample in skip mode.
    if [ ! -s coverage.bed ]; then
        if [ "${on_empty}" = "skip" ]; then
            >&2 echo "WARNING: ${empty_message}"
        else
            >&2 echo "ERROR: ${empty_message}"
            >&2 echo "Set coverage.on_empty: skip to drop this sample."
            exit 1
        fi
    fi
    """
}
