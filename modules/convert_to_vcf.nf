process CONVERT_TO_VCF {
    cpus params.consensus.cpu
    maxForks 1
    afterScript 'stage_cleanup.sh'

    publishDir "${params.outdir}/", mode: 'move', pattern: 'consensus.vcf.gz', enabled: params.output.format == 'vcf' && params.output.type == 'single'

    input:
    path('tmp/consensus.bcf')

    output:
    path("consensus.vcf.gz"), emit: consensus

    script:
    """
    bcftools view --threads ${task.cpus} -Oz9 -o consensus.vcf.gz tmp/consensus.bcf
    """

    stub:
    """
    touch consensus.vcf.gz
    """
}