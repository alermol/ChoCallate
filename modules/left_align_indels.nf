process LEFT_ALIGN_INDELS {
    maxForks params.bam_preparation.left_align_indels.forks
    cpus params.bam_preparation.left_align_indels.cpu
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
    gatk LeftAlignIndels -I tmp/input.bam -O output.bam -R tmp/ref_genome.fasta --sequence-dictionary tmp/ref_genome.dict -OBI false ${params.bam_preparation.left_align_indels.extra_args}
    """

    stub:
    """
    touch output.bam
    """
}