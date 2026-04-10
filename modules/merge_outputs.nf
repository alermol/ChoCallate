process MERGE_OUTPUTS {
    cpus params.consensus.cpu
    maxForks 1
    afterScript 'stage_cleanup.sh'
    
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.bcf', enabled: params.output.type == 'single' && params.output.format == 'bcf'

    input:
    path('tmp/?.bcf', arity: '1..*')

    output:
    path("consensus.bcf"), emit: consensus

    script:
    def regex = "ALT ~ \"\\.\""
    """
    for file in tmp/*.bcf; do bcftools index --threads ${task.cpus} --csi \$file; done

    bcftools merge --force-single --force-samples --threads ${task.cpus} -Ou tmp/*.bcf | bcftools norm -m -any --threads ${task.cpus} -Ou | bcftools filter --threads ${task.cpus} -e '${regex}' -Ou 2>/dev/null | bcftools norm --threads ${task.cpus} -m +any -Ob -o consensus.bcf
    """
}