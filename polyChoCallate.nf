// ChoCallate
// Pipeline for SNV/INDEL calling using *Cho*rus of *Call*ers

params.samples_tsv = 'input.tsv'
params.outdir = 'ChoCallate_output'
params.min_coverage = 5
params.min_base_quality = 20
params.samtools_min_map_qual = 10
params.min_snp_qual = 20
params.ploidy = 2
params.reads_type = 'pe' // se - single-end reads; pe - pair-end reads
params.reads_source = 'gbs' // gbs - Genotyping-by-sequencing; wgs - Whole Genome Sequencing
params.bowtie2_cpu = 10
params.bowtie2_forks = 1
params.freebayes_forks = 1
params.gatk4_forks = 1
params.snver_forks = 1


// Main workflow definition
workflow {
    // Create channel from samples TSV file and parse into tuples
    Channel
        .fromPath( params.samples_tsv )
        .splitCsv( header: false, sep: '\t' )
        .map { row -> tuple( row[0], file(row[1]), file(row[2]) ) }
        .set { sample_run_ch }

    // Define reference files
    ref_index = file(params.reference_index)
    ref_genome = file(params.reference_genome)

    // Create reference genome index files
    create_faidx(ref_genome)
    create_sequence_dictionary(ref_genome)

    // Map reads to reference genome
    map_reads(sample_run_ch, ref_index)

    // Left-align indels in BAM files
    left_align_indels(map_reads.out.bam, 
                      create_sequence_dictionary.out.gen_dict, 
                      create_faidx.out.ref_genome)

    // Index the aligned BAM files
    bam_indexing(left_align_indels.out.lai_bam)

    // Generate coverage information from BAM files
    bam_cov_generation(bam_indexing.out.ind_bam.map{it[0]})

    // Call SNVs/INDELs using different callers
    freebayes_snps = freebayes_calling(bam_indexing.out.ind_bam, 
                                      bam_cov_generation.out.coverage, 
                                      create_faidx.out.ref_genome)
    gatk4_snps = gatk4_calling(bam_indexing.out.ind_bam, 
                              bam_cov_generation.out.coverage, 
                              create_faidx.out.ref_genome,
                              create_sequence_dictionary.out.gen_dict)
    snver_snps = snver_calling(bam_indexing.out.ind_bam, 
                              bam_cov_generation.out.coverage, 
                              create_faidx.out.ref_genome)

    // Process SNPs
    merged_snps_vcfs = gatk4_snps.snps_vcf
        .join(freebayes_snps.snps_vcf)
        .join(snver_snps.snps_vcf)
        .map { sample, gatk, freebayes, snver ->
            tuple(sample, gatk, freebayes, snver)
        }
    generate_final_vcf_snps(merged_snps_vcfs)
    process_final_vcf_snps(generate_final_vcf_snps.out.fvcf, create_faidx.out.ref_genome.map{it[1]})

    // Process INDELs
    merged_indels_vcfs = gatk4_snps.indels_vcf
        .join(freebayes_snps.indels_vcf)
        .join(snver_snps.indels_vcf)
        .map { sample, gatk, freebayes, snver ->
            tuple(sample, gatk, freebayes, snver)
        }
    generate_final_vcf_indels(merged_indels_vcfs)
    process_final_vcf_indels(generate_final_vcf_indels.out.fvcf, create_faidx.out.ref_genome.map{it[1]})

    if (!params.debug) {
        cleanup(map_reads.out.bam,
                bam_indexing.out.ind_bam,
                bam_cov_generation.out.coverage,
                generate_final_vcf_snps.out.fvcf,
                merged_snps_vcfs,
                process_final_vcf_snps.out.final_vcf_snp,
                process_final_vcf_indels.out.final_vcf_indel)
    }
}

// Cleanup temporary files after workflow completion
workflow.onComplete {
    def workDir = workflow.workDir ? file(workflow.workDir) : null
    if (workDir?.exists() && !params.debug) {
        workDir.deleteDir()
        }
}

// Process to create FASTA index file (for FreeBayes)
process create_faidx {
    maxForks 1
    cpus 1
    
    tag "${ref_genome}-faidx"

    input:
    path(ref_genome)

    output:
    tuple path("${ref_genome}"), path("${ref_genome}.fai"), emit: ref_genome

    script:
    """
    samtools faidx --threads ${task.cpus} ${ref_genome}
    """
}

// Process to create sequence dictionary (for GATK4)
process create_sequence_dictionary {
    maxForks 1

    tag "${ref_genome.baseName}-CreateSequenceDictionary"

    input:
    path(ref_genome)
    
    output:
    path("${ref_genome.baseName}.dict"), emit: gen_dict
    
    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    """
}

// Process to map reads to reference genome
process map_reads {
    maxForks params.bowtie2_forks
    cpus params.bowtie2_cpu

    tag "${sample_id}-bowtie2"
    
    input:
    tuple val(sample_id), val(read1), val(read2)
    val(ref_index)

    output:
    path("${sample_id}.tmp_bam"), emit: bam

    script:
    if ( params.reads_type == 'pe' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -1 ${read1} -2 ${read2} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.samtools_min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.tmp_bam
        """
    else if ( params.reads_type == 'se' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -U ${read1} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.samtools_min_map_qual} -F 4 - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.tmp_bam
        """
    else
        error 'Invalid reads type: ${params.reads_type}. Available types: se, pe'
}

// Process to left-align indels in BAM files
process left_align_indels {
    cpus 1
    maxForks 1

    tag "${bam.baseName}-LeftAlignIndels"

    input:
    path(bam)
    path(genome_dictionary)
    tuple path(ref_genome), path(genome_fai)
    
    output:
    path("${bam.baseName}.bam"), emit: lai_bam
    
    script:
    """
    gatk LeftAlignIndels -I ${bam} -O ${bam.baseName}.bam -R ${ref_genome} -OBI false
    """
}

// Process to index BAM files
process bam_indexing {
    maxForks 1
    cpus 1

    tag "${bam.baseName}-samtools"
    
    input:
    path(bam)
    
    output:
    tuple path("${bam}"), path("${bam}.csi"), emit: ind_bam
    
    script:
    """
    samtools index --csi --threads ${task.cpus} ${bam.baseName}.bam
    """
}

// Process to generate coverage information
process bam_cov_generation {
    cpus 1
    maxForks 1
    
    tag "${bam.baseName}-coverage"

    input:
    path(bam)

    output:
    path("${bam.baseName}.bed"), emit: coverage

    script:
    """
    samtools depth --threads ${task.cpus} ${bam} | \
        awk '\$3 >= ${params.min_coverage} {print \$1,\$2-1,\$2}' | \
        bedops --merge - > ${bam.baseName}.bed
    """
}

// Process to call variants using FreeBayes
process freebayes_calling {
    maxForks params.freebayes_forks

    tag "${bam.baseName}-freebayes"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_freebayes"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_freebayes"), emit: indels_vcf
    
    script:
    if ( params.reads_source == 'gbs' )
        """
        freebayes --fasta-reference ${ref_genome} --targets ${coverage} --dont-left-align-indels --ploidy ${params.ploidy} \
            --use-best-n-alleles 4 --min-alternate-qsum ${params.min_base_quality} --hwe-priors-off --no-population-priors \
            --binomial-obs-priors-off --allele-balance-priors-off --min-base-quality ${params.min_base_quality} \
            --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam ${bam} --limit-coverage 250 | \
            bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
            bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
            bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.freebayes

        bcftools view -v snps -Ov -o ${bam.baseName}.snps_freebayes ${bam.baseName}.freebayes
        
        bcftools view -v indels -Ov -o ${bam.baseName}.indels_freebayes ${bam.baseName}.freebayes
        """
    
    else if ( params.reads_source == 'wgs' )
        """
        freebayes --fasta-reference ${ref_genome} --targets ${coverage} --dont-left-align-indels --ploidy ${params.ploidy} \
            --use-best-n-alleles 4 --min-alternate-qsum ${params.min_base_quality} --hwe-priors-off --no-population-priors \
            --allele-balance-priors-off --min-base-quality ${params.min_base_quality} \
            --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam ${bam} --limit-coverage 250 | \
            bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
            bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
            bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.freebayes

        bcftools view -v snps -Ov -o ${bam.baseName}.snps_freebayes ${bam.baseName}.freebayes
        
        bcftools view -v indels -Ov -o ${bam.baseName}.indels_freebayes ${bam.baseName}.freebayes
        """
    
    else
        error 'Invalid reads source: ${params.reads_source}. Available sources: gbs, wgs'
}


// Process to call variants using GATK4
process gatk4_calling {
    maxForks params.gatk4_forks

    tag "${bam.baseName}-gatk4"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    path(ref_genome_dict)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_gatk"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_gatk"), emit: indels_vcf
    
    script:
    """
    gatk HaplotypeCaller -I ${bam} -R ${ref_genome} -mbq ${params.min_base_quality} -O ${bam.baseName}.gatk1 -L ${coverage} -ERC BP_RESOLUTION -ploidy ${params.ploidy}
    bcftools filter ${bam.baseName}.gatk1 -e'QUAL<${params.min_snp_qual}' | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools sort - | \
        bcftools view -AA --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.gatk

    bcftools view -V indels,mnps,bnd,other -Ov -o ${bam.baseName}.snps_gatk ${bam.baseName}.gatk
    
    bcftools view -v indels -Ov -o ${bam.baseName}.indels_gatk ${bam.baseName}.gatk
    """
}


// Process to call variants using SNVer
process snver_calling {
    maxForks params.snver_forks

    tag "${bam.baseName}-snver"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_snver"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_snver"), emit: indels_vcf
    
    script:
    """
    ln -s ${ref_genome} reference.fasta
    
    samtools faidx reference.fasta
    
    snver -i ${bam} -r reference.fasta -o ${bam.baseName} -l ${coverage} -bq ${params.min_base_quality} -n ${params.ploidy}
    
    bcftools reheader -f reference.fasta.fai ${bam.baseName}.filter.vcf | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.snps_snver
        
    bcftools reheader -f reference.fasta.fai ${bam.baseName}.indel.filter.vcf | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.indels_snver
    """
}


// Process to generate final consensus VCF from all callers using majority rule
// Majority rule - variant is true if detected more than 2 callers
process generate_final_vcf_indels {
    maxForks 1

    tag "${sample}-generate"

    input:
    tuple val(sample), path(vcf1), path(vcf2), path(vcf3)

    output:
    path("${sample}.vcf"), emit: fvcf

    script:
    """
    python3 ${projectDir}/scripts/process_indels_poly.py --vcf1 "${vcf1}" --vcf2 "${vcf2}" --vcf3 "${vcf3}" --sample "${sample}"
    """
}

process generate_final_vcf_snps {
    maxForks 1

    tag "${sample}-generate"

    input:
    tuple val(sample), path(vcf1), path(vcf2), path(vcf3)

    output:
    path("${sample}.vcf"), emit: fvcf

    script:
    """
    python3 ${projectDir}/scripts/process_snps_poly.py --vcf1 "${vcf1}" --vcf2 "${vcf2}" --vcf3 "${vcf3}" --sample "${sample}"
    """
}

// Process to finalize and compress the VCF file
process process_final_vcf_snps {
    maxForks 1
    cpus 1

    tag "${vcf.baseName}-finalizeSNP"
    publishDir "${params.outdir}/${vcf.baseName}/", mode: 'copy', pattern: '*.snps.vcf.gz'

    input:
    path(vcf)
    path(ref_genome_fai)

    output:
    path("${vcf.baseName}.snps.vcf.gz"), emit: final_vcf_snp

    script:
    """
    bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} ${vcf} | bcftools sort -Oz -o ${vcf.baseName}.snps.vcf.gz
    """
}

process process_final_vcf_indels {
    maxForks 1
    cpus 1

    tag "${vcf.baseName}-finalizeINDEL"
    publishDir "${params.outdir}/${vcf.baseName}/", mode: 'copy', pattern: '*.indels.vcf.gz'

    input:
    path(vcf)
    path(ref_genome_fai)

    output:
    path("${vcf.baseName}.indels.vcf.gz"), emit: final_vcf_indel

    script:
    """
    bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} ${vcf} | bcftools sort -Oz -o ${vcf.baseName}.indels.vcf.gz
    """
}

process cleanup {
    maxForks 1

    tag "${bam.baseName}-cleanup"

    input:
    path(tmp_bam)
    tuple path(bam), path(index)
    path(coverage)
    path(snps_vcf)
    tuple val(sample), path(snp_vcf1), path(snp_vcf2), path(snp_vcf3)
    path(final_vcf_snps)
    path(final_vcf_indels)


    script:
    """
    for i in ${tmp_bam} ${bam} ${index} ${coverage} ${snp_vcf1} ${snp_vcf2} \
        ${snp_vcf3} ${snps_vcf}; do rm -r \$(dirname \$(realpath \${i})); done
    """
}

