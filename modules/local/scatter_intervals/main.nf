process SCATTER_INTERVALS {
    cpus params.calling.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

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
    """
    mkdir -p tmp/bed_chunks chunks

    gatk SplitIntervals \
        --java-options "-Djava.io.tmpdir=\$TMPDIR" \
        --tmp-dir \$TMPDIR \
        -R tmp/ref_genome.fasta \
        -L tmp/coverage.bed \
        --scatter-count ${scatter_count} \
        -O tmp/bed_chunks/

    parallel -j ${task.cpus} \
        'gatk IntervalListToBed -I {} -O {//}/{/.}.bed.tmp
         cut -f 1-3 {//}/{/.}.bed.tmp > chunks/{/.}.bed
         rm -f {} {//}/{/.}.bed.tmp' \
        ::: tmp/bed_chunks/*
    """
}
