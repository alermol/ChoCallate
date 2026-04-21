process MERGE_OUTPUTS {
    cpus params.consensus.cpu
    maxForks 1
    afterScript 'stage_cleanup.sh'
    
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.bcf', enabled: params.output.type == 'single' && params.output.format == 'bcf'

    input:
    path('tmp/?.bcf', arity: '1..*')
    path('tmp/coverage.bed')

    output:
    path("consensus.bcf"), emit: consensus

    script:
    def regex = 'ALT~\\"\\.\\"'
    """
    parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf

    mkdir -p tmp/bed_chunks/
    mkdir -p tmp/vcf_chunks/
    bedops --chop 10000 tmp/coverage.bed > tmp/coverage.bed.chopped
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed.chopped tmp/bed_chunks/

    parallel -j ${task.cpus} \
    "bcftools merge --force-single --force-samples --regions-file {} --threads 1 -Ou tmp/*.bcf | bcftools norm -m -any --threads 1 -Ou | bcftools filter --threads 1 -e '${regex}' -Ou | bcftools norm --threads 1 -m +any -Ob -o tmp/vcf_chunks/{#}.bcf" ::: tmp/bed_chunks/*.bed

    bcftools concat --naive --threads ${task.cpus} -Ob tmp/vcf_chunks/*.bcf | bcftools sort -Ob -o "consensus.bcf"
    """
}