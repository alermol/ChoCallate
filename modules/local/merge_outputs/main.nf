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
    // Monomorphic records must survive the merge when `remove_invariant` is not
    // set: per-sample consensus files contain all sites (see the all-sites
    // invariant) and filtering ALT="." here dropped exactly those sites.
    // A pass-through `bcftools view -Ou` keeps the pipeline shape identical.
    def invariant_filter = params.output.remove_invariant ? "bcftools filter --threads ${task.cpus} -e 'ALT~\"\\.\" || COUNT(GT=\"RR\")=N_SAMPLES' -Ou" : "bcftools view --threads ${task.cpus} -Ou"
    def output_format = params.output.format == 'vcf' ? "-Oz -o consensus.vcf.gz" : "-Ob -o consensus.bcf"
    def split_multiallelic = params.output.split_multiallelic ? "" : "| bcftools norm --threads ${task.cpus} -m +any -Ou"
    def fill_tags = "AN,AC,AF,NS,AC_Hom,AC_Het,MAF,TYPE,F_MISSING,'DP:1=int(sum(FORMAT/DP))'"
    """
    n_bcf=\$(ls -1 tmp/*.bcf | wc -l)

    if [ "\$n_bcf" -gt 1 ]; then
        parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf
    fi

    if [ "\$n_bcf" -eq 1 ]; then
        ${invariant_filter} tmp/1.bcf \
        | bcftools +fill-tags -Ou -- -t ${fill_tags} \
        | bcftools view ${output_format}
    else
        ls -1v tmp/*.bcf > tmp/merge_list.txt
        merge_bcf_tree.py --cpus ${task.cpus} --file-list tmp/merge_list.txt \
        | bcftools norm -m -any --threads ${task.cpus} -Ou \
        | drop_refonly_duplicates.py \
        | ${invariant_filter} ${split_multiallelic} \
        | fill_missing_ad.py \
        | bcftools +fill-tags -Ou -- -t ${fill_tags} \
        | bcftools sort ${output_format}
    fi
    """
}
