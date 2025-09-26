process CLEANUP_SAMPLE_TEMP {
    tag "${sample}"
    
    input:
    tuple val(sample), path(output_vcf_snps)
    tuple val(sample), path(output_vcf_indels)
    tuple path(bam), path(bam_index)
    path(coverage_bed)
    path(zero_bcf)
    tuple val(sample), path("${sample}.snps_*.bcf", arity: '3..*')
    tuple val(sample), path("${sample}.indels_*.bcf", arity: '3..*')
    val(merged)
    val(cleanup_intermediate_bam)
    val(cleanup_intermediate_bcf)
    
    script:
    """
    log_msg() {
        local level="\$1"
        local msg="\$2"
        echo "[\$(date -Iseconds)] [${sample}] [CLEANUP_SAMPLE_TEMP] [\${level}] \${msg}"
    }

    log_msg "INFO" "Starting sample cleanup process (symlink-aware)"

    remove_file_follow_symlink() {
        local f="\$1"
        if [ -L "\${f}" ]; then
            local target
            target=\$(realpath "\${f}")
            if [ -e "\$target" ]; then
                rm -f "\$target" && log_msg "INFO" "Cleaned up symlink target: \$target"
            fi
            rm -f "\$f" && log_msg "INFO" "Cleaned up symlink: \$f"
        elif [ -e "\$f" ]; then
            rm -f "\$f" && log_msg "INFO" "Cleaned up file: \$f"
        fi
    }

    for vcf in ${output_vcf_snps}; do
        [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        [ -e "\${vcf}.csi" ] && remove_file_follow_symlink "\${vcf}.csi"
    done

    for vcf in ${output_vcf_indels}; do
        [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        [ -e "\${vcf}.csi" ] && remove_file_follow_symlink "\${vcf}.csi"
    done

    if [ "${cleanup_intermediate_bam}" = "true" ]; then
        [ -e "${bam}" ] && remove_file_follow_symlink "${bam}"
        [ -e "${bam_index}" ] && remove_file_follow_symlink "${bam_index}"
        log_msg "INFO" "Cleaned up BAM files"
    fi

    [ -e "${coverage_bed}" ] && remove_file_follow_symlink "${coverage_bed}"
    [ -e "${zero_bcf}" ] && remove_file_follow_symlink "${zero_bcf}"

    if [ "${cleanup_intermediate_bcf}" = "true" ]; then
        for vcf in "${sample}.snps_"*.bcf "${sample}.indels_"*.bcf; do
            [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        done
        log_msg "INFO" "Cleaned up consensus VCFs"
    fi

    for symlink in ./*; do
        if [ -L "\${symlink}" ]; then
            rm -f "\${symlink}" && log_msg "INFO" "Cleaned up symlink in current directory: \${symlink}"
        fi
    done

    log_msg "INFO" "Sample cleanup completed"
    """
}


