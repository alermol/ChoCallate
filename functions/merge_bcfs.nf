process MERGE_BCFS {
    maxForks params.merge_bcfs_forks
    cpus params.merge_bcfs_cpus

    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.snps.bcf', enabled: !params.output_vcf
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.indels.bcf', enabled: !params.output_vcf
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.snps.vcf.gz', enabled: params.output_vcf
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.indels.vcf.gz', enabled: params.output_vcf

    input:
    path('?.snps', arity: '1..*')
    path('?.indels', arity: '1..*')

    output:
    path 'final.snps.*', emit: final_snps
    path 'final.indels.*', emit: final_indels
    val true, emit: merged

    script:
    """
    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Process started - Merging BCFs"

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Indexing SNPs BCFs"
    parallel -j ${task.cpus} 'bcftools index --csi --threads 1 {}' ::: *.snps
    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Indexing indels BCFs"
    parallel -j ${task.cpus} 'bcftools index --csi --threads 1 {}' ::: *.indels

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Merging SNPs BCFs"
    bcftools merge --force-single --threads ${task.cpus} -Ob -o final.snps.bcf ${params.bcftools_merge_extra_args} *.snps
    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Merging indels BCFs"
    bcftools merge --force-single --threads ${task.cpus} -Ob -o final.indels.bcf ${params.bcftools_merge_extra_args} *.indels

    if [ "${params.output_vcf}" = "true" ]; then
        bcftools view -Oz -o final.snps.vcf.gz final.snps.bcf
        bcftools view -Oz -o final.indels.vcf.gz final.indels.bcf
        rm -f final.snps.bcf final.indels.bcf
    fi

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Cleaning up intermediate BCFs"
    find . -name '*.snps' -type l -delete
    find . -name '*.indels' -type l -delete
    find . -name '*.csi' -type f -delete

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Process completed - Merged BCFs"
    """
}


