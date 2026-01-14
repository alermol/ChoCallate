process GENERATE_CONSENSUS {
    maxForks params.cons_forks
    cpus params.cons_cpus

    tag "${sample}"

    publishDir "${params.outdir}/per_sample/${sample}/", mode: 'copy', pattern: '*.snps.bcf', enabled: !params.output_vcf && params.per_sample_out
    publishDir "${params.outdir}/per_sample/${sample}/", mode: 'copy', pattern: '*.indels.bcf', enabled: !params.output_vcf && params.per_sample_out
    publishDir "${params.outdir}/per_sample/${sample}/", mode: 'copy', pattern: '*.snps.vcf.gz', enabled: params.output_vcf && params.per_sample_out
    publishDir "${params.outdir}/per_sample/${sample}/", mode: 'copy', pattern: '*.indels.vcf.gz', enabled: params.output_vcf && params.per_sample_out

    input:
    tuple val(sample), path("${sample}.snps_*.bcf", arity: '3..*')
    tuple val(sample), path("${sample}.indels_*.bcf", arity: '3..*')
    path(ref_genome_fai)
    path(zero_bcf)
    val(cons_threshold)

    output:
    tuple val("${sample}"), path("${sample}.snps.*"), emit: final_snps
    tuple val("${sample}"), path("${sample}.indels.*"), emit: final_indels

    script:
    """
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Process started - Generating consensus BCFs"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Consensus threshold: ${cons_threshold}"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Window size: ${params.win_size}"
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Creating intermediate subfolder for consensus generation"
    mkdir -p "${sample}_consensus_gen"
    pushd "${sample}_consensus_gen" > /dev/null
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Creating symlinks to input files"
    ln -sf "../${ref_genome_fai}" ref_genome.fasta.fai
    ln -sf "../${zero_bcf}" zero.bcf
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Creating symlinks to BCF files"
    for bcf in ../${sample}.snps_*.bcf; do
        ln -sf "\$bcf" "\$(basename \$bcf)"
    done
    
    for bcf in ../${sample}.indels_*.bcf; do
        ln -sf "\$bcf" "\$(basename \$bcf)"
    done
    
    awk -v OFS='\t' '{print \$1,"0",\$2}' ref_genome.fasta.fai > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Processing SNPs..."
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indexing SNP BCF files with bcftools"
    for i in ${sample}.snps_*; do bcftools index \${i}; done

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indexing zero BCF with bcftools"
    bcftools index zero.bcf

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Extracting chromosome names from zero BCF"
    bcftools query -f '%CHROM' zero.bcf | uniq > zero_chr_names.txt

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Running consensus generation for SNPs in parallel (${task.cpus} threads)"
    parallel -j ${task.cpus} 'consensus_generation.sh {1} {#} ${sample} "snps" ${cons_threshold}' :::: genome_intervals.bed
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Collecting generated BCF files for concatenation"
    find all_chrs/ -name '*.bcf' -type f > vcf_files.txt
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Concatenating, reheading, and sorting SNP consensus BCFs"
    bcftools concat --naive -Ob --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ref_genome.fasta.fai | \
        bcftools sort -Ob -o ${sample}.snps.bcf
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] SNPs consensus BCF created"
    rm -r all_chrs/*

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Processing indels..."
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indexing indel BCF files with bcftools"
    for i in ${sample}.indels_*; do bcftools index \${i}; done

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Running consensus generation for indels in parallel (${task.cpus} threads)"
    parallel -j ${task.cpus} 'consensus_generation.sh {1} {#} ${sample} "indels" ${cons_threshold}' :::: genome_intervals.bed

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Collecting generated indel BCF files for concatenation"
    find all_chrs/ -name '*.bcf' -type f > vcf_files.txt

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Concatenating, reheading, and sorting indel consensus BCFs"
    bcftools concat --naive -Ob --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ref_genome.fasta.fai | \
        bcftools sort -Ob -o ${sample}.indels.bcf
        
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indels consensus BCF created"
    
    if [ "${params.output_vcf}" = "false" ]; then
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Moving final outputs to parent directory"
        mv ${sample}.snps.bcf ../${sample}.snps.bcf
        mv ${sample}.indels.bcf ../${sample}.indels.bcf
    else
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Converting consensus BCF to VCF"
        bcftools view -Oz -o ../${sample}.snps.vcf.gz ${sample}.snps.bcf
        bcftools view -Oz -o ../${sample}.indels.vcf.gz ${sample}.indels.bcf
    fi
    
    popd > /dev/null

    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${sample}_consensus_gen"
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${ref_genome_fai}" "${zero_bcf}" "${sample}.snps_"*.bcf "${sample}.indels_"*.bcf 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Process completed - Consensus BCFs generated"
    """
}


