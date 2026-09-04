process MERGE_OUTPUTS {
    maxForks 1
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'
    
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.bcf', enabled: params.output.type == 'single' && params.output.format == 'bcf'
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.vcf.gz', enabled: params.output.type == 'single' && params.output.format == 'vcf'

    input:
    path('tmp/?.bcf', arity: '1..*')

    output:
    path("consensus.bcf"), optional: true
    path("consensus.vcf.gz"), optional: true

    script:
    def regex = params.output.remove_invariant ? 'ALT~\"\\.\" || COUNT(GT=\"hom\")=N_SAMPLES || COUNT(GT=\"het\")=N_SAMPLES' : 'ALT~\"\\.\"'
    def output_format = params.output.format == 'vcf' ? "-Oz -o consensus.vcf.gz" : "-Ob -o consensus.bcf"
    def split_multiallelic = params.output.split_multiallelic ? "" : "| bcftools norm --threads ${task.cpus} -m +any -Ou"
    """
    n_bcf=\$(ls -1 tmp/*.bcf | wc -l)

    if [ "\$n_bcf" -gt 1 ]; then
        parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf
    fi

    if [ "\$n_bcf" -eq 1 ]; then
        bcftools filter --threads ${task.cpus} -e '${regex}' ${output_format} tmp/1.bcf
    else
        ls -1 tmp/*.bcf > tmp/merge_list.txt
        merge_bcf_tree.py --cpus ${task.cpus} --file-list tmp/merge_list.txt \
        | bcftools norm -m -any --threads ${task.cpus} -Ou \
        | bcftools filter --threads ${task.cpus} -e '${regex}' -Ou ${split_multiallelic} \
        | fill_missing_ad.py \
        | bcftools sort ${output_format}
    fi
    """
}
