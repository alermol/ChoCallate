process MERGE_OUTPUTS {
    cpus 1
    maxForks 1
    afterScript 'stage_cleanup.sh'
    
    publishDir "${params.outdir}", mode: 'move', pattern: 'consensus.bcf', enabled: params.output.type == 'single' && params.output.format == 'bcf'

    input:
    path('tmp/?.bcf', arity: '1..*')

    output:
    path("consensus.bcf"), emit: consensus

    script:
    """
    for file in tmp/*.bcf; do bcftools index --threads ${task.cpus} --csi \$file; done

    bcftools merge --force-single --threads ${task.cpus} -Ob -o consensus.bcf tmp/*.bcf
    """

    stub:
    """
    touch consensus.bcf
    """
}