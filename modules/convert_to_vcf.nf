process CONVERT_TO_VCF_SAMPLE {
    cpus params.consensus.cpu
    maxForks 1
    //afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    publishDir "${params.output.directory}/per_sample", mode: 'move', pattern: '*.vcf.gz', enabled: params.output.format == 'vcf' && params.output.type == 'sample'

    input:
    tuple val(sample_id), path('tmp/consensus.bcf')

    output:
    path("${sample_id}.vcf.gz"), emit: consensus

    script:
    """
    bcftools view --threads ${task.cpus} -Oz9 -o ${sample_id}.vcf.gz tmp/consensus.bcf
    """
}

process CONVERT_TO_VCF_SINGLE {
    cpus params.consensus.cpu
    maxForks 1
    //afterScript 'stage_cleanup.sh'

    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.vcf.gz', enabled: params.output.format == 'vcf' && params.output.type == 'single'

    input:
    path('tmp/consensus.bcf')

    output:
    path("consensus.vcf.gz"), emit: consensus

    script:
    """
    bcftools view --threads ${task.cpus} -Oz9 -o consensus.vcf.gz tmp/consensus.bcf
    """
}