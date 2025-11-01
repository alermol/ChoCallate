process CALLING {
    maxForks params.calling_forks
    cpus 1

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    tuple path(ref_genome), path(ref_genome_fai)
    path(coverage_bed)
    val(ploidy)
    path(ref_genome_dict)
    val(effective_callers)
    val(parallel_cpus)
    val(bcftools_mpileup_extra_args)
    val(bcftools_call_extra_args)
    val(freebayes_extra_args)
    val(gatk4_extra_args)
    val(vardict_extra_args)
    val(snver_extra_args)

    output:
    tuple val("${bam.baseName}"), path("*.snps_*.bcf"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("*.indels_*.bcf"), emit: indels_vcf

    script:
    """
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Process started - Running variant callers"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Using callers: ${effective_callers}"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Parallel CPUs: ${parallel_cpus}"
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Creating intermediate subfolder for variant calling"
    mkdir -p "${bam.baseName}_calling"
    cd "${bam.baseName}_calling"
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Creating symlinks to input files"
    ln -sf "../${bam}" input.bam
    ln -sf "../${bam_index}" input.bam.csi
    ln -sf "../${ref_genome}" reference.fasta
    ln -sf "../${ref_genome_fai}" reference.fasta.fai
    ln -sf "../${coverage_bed}" coverage.bed
    ln -sf "../${ref_genome_dict}" reference.dict
    
    touch callers_commands.sh

    if [[ ",${effective_callers}," == *"bcftools"* ]]; then
        echo "bash ${projectDir}/bin/bcftools_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.bcftools_cpu} '${bcftools_mpileup_extra_args}' '${bcftools_call_extra_args}'" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"freebayes"* ]]; then
        echo "bash ${projectDir}/bin/freebayes_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.reads_source} ${params.min_base_quality} ${params.min_snp_qual} '${freebayes_extra_args}'" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"gatk"* ]]; then
        echo "bash ${projectDir}/bin/gatk4_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} '${gatk4_extra_args}'" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"vardict"* ]]; then
        echo "bash ${projectDir}/bin/vardict_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.vardict_cpu} '${vardict_extra_args}'" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"snver"* ]]; then
        echo "bash ${projectDir}/bin/snver_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} '${snver_extra_args}'" >> callers_commands.sh
    fi

    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Executing ${parallel_cpus} callers in parallel"
    parallel -j ${parallel_cpus} '{}' :::: callers_commands.sh
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Moving final outputs to parent directory"
    mv *.snps_*.bcf ../
    mv *.indels_*.bcf ../
    
    cd ..
    
    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${bam.baseName}_calling"
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${bam}" "${bam_index}" "${ref_genome}" "${ref_genome_fai}" "${coverage_bed}" "${ref_genome_dict}" 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Process completed - Variant calling finished"
    """
}


