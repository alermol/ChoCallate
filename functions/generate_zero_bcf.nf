process GENERATE_ZERO_BCF {
    maxForks params.zero_bcf_forks
    cpus params.zero_bcf_cpu

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    tuple path(ref_genome), path(ref_genome_fai)
    path(coverage_bed)

    output:
    path("zero.bcf"), emit: zero_bcf

    script:
    String genotype = (["0"] * params.ploidy).join("/")
    
    """
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Process started - Generating zero BCF"
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Creating intermediate subfolder for zero BCF generation"
    mkdir -p "${bam.baseName}_zero_bcf_gen"
    cd "${bam.baseName}_zero_bcf_gen"
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Creating symlinks to input files"
    ln -sf "../${bam}" input.bam
    ln -sf "../${bam_index}" input.bam.csi
    ln -sf "../${ref_genome}" ref_genome.fasta
    ln -sf "../${ref_genome_fai}" ref_genome.fasta.fai
    ln -sf "../${coverage_bed}" coverage.bed

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Writing minimal VCF header for bcftools reheader fallback"
    cat > minimal_header.txt <<EOF
    ##fileformat=VCFv4.3
    ##source=ChoCallateZeroVCF
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
    EOF
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Running bcftools mpileup with min base quality ${params.min_base_quality}"
    bcftools mpileup -Ov --count-orphans --fasta-ref ref_genome.fasta --threads ${task.cpus} --max-depth 1 \
        --min-BQ ${params.min_base_quality} --regions-file coverage.bed input.bam | \
        awk -v OFS='\t' -v gen=${genotype} '{if(\$0 !~ /#/) print \$1,\$2,\$3,\$4,".","100",".",".","GT",gen; else print \$0}' | \
        awk -v OFS='\t' '{if(length(\$4) == 1 || \$0 ~ /#/) print \$0}' | \
        bcftools reheader -h minimal_header.txt -f ref_genome.fasta.fai - | \
        bcftools view -Ob -o zero.bcf

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Moving final output to parent directory"
    mv zero.bcf ../zero.bcf
    
    cd ..
    
    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${bam.baseName}_zero_bcf_gen"
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${bam}" "${bam_index}" "${ref_genome}" "${ref_genome_fai}" "${coverage_bed}" 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Process completed - Zero BCF created"
    """
}


