process MERGE_BCFS {
    maxForks params.merge_bcfs_forks
    cpus params.merge_bcfs_cpus

    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.snps.bcf', enabled: !params.output_vcf && !params.merge_variants
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.indels.bcf', enabled: !params.output_vcf && !params.merge_variants

    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.snps.vcf.gz', enabled: params.output_vcf && !params.merge_variants
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.indels.vcf.gz', enabled: params.output_vcf && !params.merge_variants

    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.merged.bcf', enabled: !params.output_vcf && params.merge_variants
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.merged.vcf.gz', enabled: params.output_vcf && params.merge_variants


    input:
    path('?.snps', arity: '1..*')
    path('?.indels', arity: '1..*')
    path('?.merged', arity: '1..*')

    output:
    path 'final.snps.*'
    path 'final.indels.*'
    path 'final.merged.*'
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
        echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Converting BCF to VCF"
        bcftools view --threads ${task.cpus} -Oz -o final.snps.vcf.gz final.snps.bcf
        bcftools view --threads ${task.cpus} -Oz -o final.indels.vcf.gz final.indels.bcf
        rm -f final.snps.bcf final.indels.bcf

        echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Indexing VCF with tabix"
        tabix -C --threads ${task.cpus} final.snps.vcf.gz
        tabix -C --threads ${task.cpus} final.indels.vcf.gz
        
        echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Concatenating VCFs"
        bcftools concat -a --threads ${task.cpus} -Oz -o final.merged.vcf.gz final.snps.vcf.gz final.indels.vcf.gz
    else
        echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Indexing BCFs with tabix"
        tabix -C --threads ${task.cpus} final.snps.bcf
        tabix -C --threads ${task.cpus} final.indels.bcf

        echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Concatenating BCFs"
        bcftools concat -a --threads ${task.cpus} -Ob -o final.merged.bcf final.snps.bcf final.indels.bcf
    fi


    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Cleaning up intermediate BCFs"
    find . -name '*.snps' -type l -delete
    find . -name '*.indels' -type l -delete
    find . -name '*.csi' -type f -delete

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Process completed - Merged BCFs"
    """
}


