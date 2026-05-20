process MERGE_OUTPUTS {
    maxForks 1
    cpus { Runtime.runtime.availableProcessors() - task.attempt }
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'
    errorStrategy 'retry'
    maxRetries Runtime.runtime.availableProcessors() - 1
    
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.bcf', enabled: params.output.type == 'single' && params.output.format == 'bcf'
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.vcf.gz', enabled: params.output.type == 'single' && params.output.format == 'vcf'

    input:
    path('tmp/?.bcf', arity: '1..*')
    path('tmp/coverage.bed')

    output:
    path("consensus.bcf"), optional: true
    path("consensus.vcf.gz"), optional: true

    script:
    def regex = 'ALT~\\"\\.\\"'
    def output_format = params.output.format == 'vcf' ? "-Oz -o consensus.vcf.gz" : "-Ob -o consensus.bcf"
    """
    parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf

    mkdir -p tmp/bed_chunks/
    mkdir -p tmp/vcf_chunks/
    bedops --chop 10000 tmp/coverage.bed > tmp/coverage.bed.chopped
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed.chopped tmp/bed_chunks/

    parallel -j ${task.cpus} \
    "bcftools merge --force-single --force-samples --regions-file {} --threads 1 -Ou tmp/*.bcf | bcftools norm -m -any --threads 1 -Ou | bcftools filter --threads 1 -e '${regex}' -Ou | bcftools norm --threads 1 -m +any -Ob -o tmp/vcf_chunks/{#}.bcf" ::: tmp/bed_chunks/*.bed

    bcftools concat --naive --threads ${task.cpus} -Ob tmp/vcf_chunks/*.bcf | bcftools sort ${output_format}
    """
}
