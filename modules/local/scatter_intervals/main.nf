process SCATTER_INTERVALS {
    label 'process_medium'
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/coverage.bed")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/ref_genome.dict")

    output:
    tuple val(sample_id), path("chunks/*.bed", arity: '1..*'),
        emit: chunks

    script:
    def scatter_count = Math.max(params.calling.cpu ?: 1, params.consensus.cpu ?: 1)
    def split_heap_option = ResourceUtils.javaHeapOption(task.memory)
    def parallel_heap_option =
        ResourceUtils.javaHeapOption(task.memory, task.cpus)
    """
    mkdir -p tmp/bed_chunks chunks

    gatk SplitIntervals \
        --java-options "-Djava.io.tmpdir=\$TMPDIR ${split_heap_option}" \
        --tmp-dir "\$TMPDIR" \
        -R tmp/ref_genome.fasta \
        -L tmp/coverage.bed \
        --scatter-count ${scatter_count} \
        -O tmp/bed_chunks/

    cat > tmp/run_chunk.sh <<'EOF'
    set -euo pipefail
    dir=\${1%/*}
    stem=\${1##*/}
    stem=\${stem%.*}
    gatk IntervalListToBed \
        --java-options "-Djava.io.tmpdir=\$TMPDIR ${parallel_heap_option}" \
        --TMP_DIR "\$TMPDIR" \
        -I "\$1" -O "\$dir/\$stem.bed.tmp"
    cut -f 1-3 "\$dir/\$stem.bed.tmp" > "chunks/\$stem.bed"
    rm -f "\$1" "\$dir/\$stem.bed.tmp"
    EOF

    # SIGKILLed children appear as signal 9, exit 137, or GATK's 247.
    parallel --halt now,fail=1 --joblog tmp/parallel.log \
        -j ${task.cpus} bash tmp/run_chunk.sh {} \
        ::: tmp/bed_chunks/* || {
        status=\$?
        if awk 'NR > 1 && (\$7 == 137 || \$7 == 247 || \$8 == 9) {
            found = 1 } END { exit !found }' tmp/parallel.log; then
            exit 137
        fi
        exit "\${status}"
    }

    # A failed IntervalListToBed can still leave an empty file behind; fail
    # instead of emitting a region that no caller will process.
    for f in chunks/*.bed; do
        [ -s "\$f" ] || { echo "empty chunk \$f" >&2; exit 1; }
    done
    """
}
