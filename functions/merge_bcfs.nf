process MERGE_BCFS {
    maxForks params.merge_bcfs_forks
    cpus params.merge_bcfs_cpus

    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.snps.bcf'
    publishDir "${params.outdir}/", mode: 'move', pattern: 'final.indels.bcf'

    input:
    path('?.snps.bcf', arity: '1..*')
    path('?.indels.bcf', arity: '1..*')

    output:
    path 'final.snps.bcf', emit: final_snps
    path 'final.indels.bcf', emit: final_indels
    val true, emit: merged

    script:
    """
    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Process started - Merging BCFs"

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Indexing SNPs BCFs"
    parallel -j ${task.cpus} 'bcftools index --csi --threads 1 {}' ::: *.snps.bcf 
    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Indexing indels BCFs"
    parallel -j ${task.cpus} 'bcftools index --csi --threads 1 {}' ::: *.indels.bcf

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Merging SNPs BCFs"
    bcftools merge --force-single --threads ${task.cpus} -Ob -o final.snps.bcf ${params.bcftools_merge_extra_args} *.snps.bcf
    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Merging indels BCFs"
    bcftools merge --force-single --threads ${task.cpus} -Ob -o final.indels.bcf ${params.bcftools_merge_extra_args} *.indels.bcf

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Cleaning up intermediate BCFs"
    find . -name '*.snps.bcf' -type l -delete
    find . -name '*.indels.bcf' -type l -delete
    find . -name '*.csi' -type f -delete

    echo "[\$(date -Iseconds)] [INFO] [MERGE_BCFS] Process completed - Merged BCFs"
    """
}


