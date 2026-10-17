process MERGE_OUTPUTS {
    maxForks 1
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.bcf', enabled: params.output.type == 'single' && params.output.format == 'bcf'
    publishDir "${params.output.directory}", mode: 'move', pattern: 'consensus.vcf.gz', enabled: params.output.type == 'single' && params.output.format == 'vcf'

    input:
    // The sample ids are carried explicitly so the merge order follows
    // samples_tsv even though `dir/*` staging preserves the file names.
    tuple val(sample_ids), path('tmp/bcfs/*', arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")

    output:
    path("consensus.bcf"), optional: true
    path("consensus.vcf.gz"), optional: true

    script:
    def sample_id_list = sample_ids.split(',').collect { id -> id.trim() }
    def merge_paths = sample_id_list.collect { id -> "tmp/bcfs/${id}.bcf" }.join(' ')
    // Per-sample consensus files contain a gVCF-style reference record
    // (ALT="<*>") at every covered position.  `bcftools merge --gvcf` expands
    // those records at variant positions of other samples and maps their
    // "other allele" AD onto the real ALT, so hom-reference samples report
    // honest ALT depth instead of a fabricated 0.  Monomorphic records must
    // survive when `remove_invariant` is not set; when it is set, the explicit
    // ALT clauses remove them (the all-hom-ref clause alone is not enough when
    // some samples have no coverage and are ./. after the merge).
    def invariant_filter = params.output.remove_invariant ? "bcftools filter --threads ${task.cpus} -e 'ALT=\".\" || ALT=\"<*>\" || COUNT(GT=\"RR\")=N_SAMPLES' -Ou" : "bcftools view --threads ${task.cpus} -Ou"
    def output_format = params.output.format == 'vcf' ? "-Oz -o consensus.vcf.gz" : "-Ob -o consensus.bcf"
    def split_multiallelic = params.output.split_multiallelic ? "" : "| bcftools norm --threads ${task.cpus} -m +any -Ou"
    def fill_tags = "AN,AC,AF,NS,AC_Hom,AC_Het,MAF,TYPE,F_MISSING,'DP:1=int(sum(FORMAT/DP))'"
    """
    n_bcf=\$(ls -1 tmp/bcfs/*.bcf | wc -l)

    if [ "\$n_bcf" -gt 1 ]; then
        parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/bcfs/*.bcf
    fi

    if [ "\$n_bcf" -eq 1 ]; then
        ${invariant_filter} tmp/bcfs/*.bcf \
        | bcftools +fill-tags -Ou -- -t ${fill_tags} \
        | bcftools view ${output_format}
    else
        printf '%s\n' ${merge_paths} > tmp/merge_list.txt
        merge_bcf_tree.py --gvcf tmp/ref_genome.fasta --cpus ${task.cpus} --file-list tmp/merge_list.txt \
        | bcftools norm -m -any --threads ${task.cpus} -Ou \
        | drop_refonly_duplicates.py \
        | ${invariant_filter} ${split_multiallelic} \
        | bcftools +fill-tags -Ou -- -t ${fill_tags} \
        | bcftools sort ${output_format}
    fi
    """
}
