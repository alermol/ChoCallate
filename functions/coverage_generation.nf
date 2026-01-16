process COVERAGE_GENERATION {
    cpus 1
    maxForks 1

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    path(custom_bed)

    output:
    path("${bam.baseName}.bed"), emit: coverage

    script:
    """
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Process started - Generating coverage information"
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Creating intermediate subfolder for coverage generation"
    mkdir -p "${bam.baseName}_coverage_gen"
    cd "${bam.baseName}_coverage_gen"
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Creating symlinks to input files"
    ln -sf "../${bam}" input.bam
    ln -sf "../${bam_index}" input.bam.csi
    mv "../${custom_bed}" custom.bed

    if [ -f "${custom_bed}" ]; then
        echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Running samtools depth with min coverage ${params.min_coverage} and custom bed file"
        samtools depth -J --threads ${task.cpus} input.bam | \
            awk '\$3 >= ${params.min_coverage} {print \$1,\$2-1,\$2}' | \
            bedops --merge - | bedops --element-of 1 - custom.bed > ${bam.baseName}.bed
    else
        echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Running samtools depth with min coverage ${params.min_coverage}"
        samtools depth -J --threads ${task.cpus} input.bam | \
            awk '\$3 >= ${params.min_coverage} {print \$1,\$2-1,\$2}' | \
            bedops --merge - > ${bam.baseName}.bed
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Moving final output to parent directory"
    mv ${bam.baseName}.bed ../${bam.baseName}.bed
    
    cd ..
    
    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${bam.baseName}_coverage_gen"
        echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${bam}" "${bam_index}" 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Process completed - Coverage BED file created"
    """
}


