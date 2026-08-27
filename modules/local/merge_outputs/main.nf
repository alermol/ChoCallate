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
    def regex = params.output.remove_invariant ? 'ALT~\"\\.\" || N_PASS(GT=\"hom\")==N_SAMPLES' : 'ALT~\"\\.\"'
    def output_format = params.output.format == 'vcf' ? "-Oz -o consensus.vcf.gz" : "-Ob -o consensus.bcf"
    def split_multiallelic = params.output.split_multiallelic ? "" : "| bcftools norm --threads ${task.cpus} -m +any -Ou"
    """
    mkdir -p tmp/merge1/
    ls tmp/*.bcf > tmp/merge1/list_bcf.txt
    n_bcf=\$(wc -l < tmp/merge1/list_bcf.txt)

    if [ "\$n_bcf" -gt 1 ]; then
        parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf
    fi

    if [ "\$n_bcf" -eq 1 ]; then
        bcftools filter --threads ${task.cpus} -e '${regex}' ${output_format} tmp/1.bcf
    elif [ "\$n_bcf" -le ${task.cpus} ]; then
        bcftools merge --force-samples --threads ${task.cpus} -Ou tmp/*.bcf \
        | bcftools norm -m -any --threads ${task.cpus} -Ou \
        | bcftools filter --threads ${task.cpus} -e '${regex}' -Ou ${split_multiallelic} | bcftools sort ${output_format}
    else
        split -n l/${task.cpus} --additional-suffix=".txt" -a 4 -d tmp/merge1/list_bcf.txt tmp/merge1/chunk_

        parallel -j ${task.cpus} \
        "bcftools merge --force-samples --force-single --threads 1 -Ou -m all -l {} -o tmp/merge1/merged_{#}.bcf" ::: tmp/merge1/chunk_*.txt

        parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/merge1/merged_*.bcf

        bcftools merge --force-samples --threads ${task.cpus} -m all -Ou tmp/merge1/merged_*.bcf \
        | bcftools norm -m -any --threads ${task.cpus} -Ou \
        | bcftools filter --threads ${task.cpus} -e '${regex}' -Ou ${split_multiallelic} | bcftools sort ${output_format}
    fi
    """
}
