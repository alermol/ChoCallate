process SCATTER_INTERVALS {
    label 'process_medium'
    cpus params.calling.cpu
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
    def heap_option = ResourceUtils.javaHeapOption(task.memory, task.cpus)
    """
    mkdir -p tmp/bed_chunks chunks

    gatk SplitIntervals \
        --java-options "-Djava.io.tmpdir=\$TMPDIR ${heap_option}" \
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
    gatk IntervalListToBed --java-options "${heap_option}" \
        -I "\$1" -O "\$dir/\$stem.bed.tmp"
    cut -f 1-3 "\$dir/\$stem.bed.tmp" > "chunks/\$stem.bed"
    rm -f "\$1" "\$dir/\$stem.bed.tmp"
    EOF

    parallel -j ${task.cpus} bash tmp/run_chunk.sh {} \
        ::: tmp/bed_chunks/*

    # A failed IntervalListToBed can still leave an empty file behind; fail
    # instead of emitting a region that no caller will process.
    for f in chunks/*.bed; do
        [ -s "\$f" ] || { echo "empty chunk \$f" >&2; exit 1; }
    done
    """
}
