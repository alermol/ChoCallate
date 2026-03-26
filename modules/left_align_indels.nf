process LEFT_ALIGN_INDELS {
    maxForks 1
    cpus 1
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path('tmp/ref_genome.fasta')
    path('tmp/ref_genome.dict')
    path('tmp/ref_genome.fasta.fai')

    output:
    tuple val(sample_id), path("output.bam"), emit: bam

    script:
    """
    gatk LeftAlignIndels -I tmp/input.bam -O output.bam -R tmp/ref_genome.fasta --sequence-dictionary tmp/ref_genome.dict -OBI false
    """

    stub:
    """
    touch output.bam
    """
}