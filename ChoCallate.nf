#!/usr/bin/env nextflow


workflow CALLING {
    take:
    indexed_bam
    coverage_file
    ref_genome
    ref_genome_dict
    ploidy

    main:
    if (ploidy == 2) {
        FREEBAYES(indexed_bam, coverage_file, ref_genome, ploidy)
        GATK4(indexed_bam, coverage_file, ref_genome, ref_genome_dict, ploidy)
        SNVER(indexed_bam, coverage_file, ref_genome, ploidy)
        BCFTOOLS(indexed_bam, coverage_file, ref_genome, ploidy)
        VARDICT(indexed_bam, coverage_file, ref_genome, ploidy)

        merged_snps_vcfs = BCFTOOLS.out.snps_vcf
            .join(FREEBAYES.out.snps_vcf)
            .join(GATK4.out.snps_vcf)
            .join(VARDICT.out.snps_vcf)
            .join(SNVER.out.snps_vcf)
            .map { sample, bcftools, freebayes, gatk, vardict, snver ->
                tuple(sample, bcftools, freebayes, gatk, vardict, snver)
            }
        merged_indels_vcfs = BCFTOOLS.out.indels_vcf
            .join(FREEBAYES.out.indels_vcf)
            .join(GATK4.out.indels_vcf)
            .join(VARDICT.out.indels_vcf)
            .join(SNVER.out.indels_vcf)
            .map { sample, bcftools, freebayes, gatk, vardict, snver ->
                tuple(sample, bcftools, freebayes, gatk, vardict, snver)
            }
    } else {
        FREEBAYES(indexed_bam, coverage_file, ref_genome, ploidy)
        GATK4(indexed_bam, coverage_file, ref_genome, ref_genome_dict, ploidy)
        SNVER(indexed_bam, coverage_file, ref_genome, ploidy)

        merged_snps_vcfs = GATK4.out.snps_vcf
            .join(FREEBAYES.out.snps_vcf)
            .join(SNVER.out.snps_vcf)
            .map { sample, gatk, freebayes, snver ->
                tuple(sample, gatk, freebayes, snver)
            }
        merged_indels_vcfs = GATK4.out.indels_vcf
            .join(FREEBAYES.out.indels_vcf)
            .join(SNVER.out.indels_vcf)
            .map { sample, gatk, freebayes, snver ->
                tuple(sample, gatk, freebayes, snver)
            }
    }

    emit:
    merged_snps_vcfs
    merged_indels_vcfs
}


// Main workflow definition
workflow {
    // Create channel from samples TSV file and parse into tuples
    Channel
        .fromPath(params.samples_tsv)
        .splitCsv(header: false, sep: '\t')
        .map{row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3]))}
        .set{sample_run_ch}

    // Define reference files
    ref_index = file(params.reference_index)
    ref_genome = file(params.reference_genome)

    // Create reference genome index files
    CREATE_FAI_INDEX(ref_genome)
    CREATE_SEQ_DICT(ref_genome)

    // Map reads to reference genome
    BOWTIE2_MAPPING(sample_run_ch, ref_index)

    // Left-align indels in BAM files
    LEFT_ALIGN_INDELS(BOWTIE2_MAPPING.out.bam,
                      CREATE_SEQ_DICT.out.gen_dict,
                      CREATE_FAI_INDEX.out.fai_index)

    // Index the aligned BAM files
    INDEXING_BAM(LEFT_ALIGN_INDELS.out.lai_bam)

    // Generate coverage information from BAM files
    COVERAGE_GENERATION(INDEXING_BAM.out.ind_bam.map{it[0]})

    // Call SNVs/INDELs using different callers
    CALLING(INDEXING_BAM.out.ind_bam,
            COVERAGE_GENERATION.out.coverage,
            CREATE_FAI_INDEX.out.fai_index,
            CREATE_SEQ_DICT.out.gen_dict,
            params.ploidy)

    // Generate final consensus
    if (params.ploidy == 2) {
        GENERATE_CONSENSUS_DIPLOID(CALLING.out.merged_snps_vcfs,
                                   CALLING.out.merged_indels_vcfs,
                                   CREATE_FAI_INDEX.out.fai_index.map{it[1]})
    } else {
        GENERATE_CONSENSUS_POLYPLOID(CALLING.out.merged_snps_vcfs,
                                     CALLING.out.merged_indels_vcfs,
                                     CREATE_FAI_INDEX.out.fai_index.map{it[1]})
    }

}

// Cleanup temporary files after workflow completion
workflow.onComplete {
    def workDir = workflow.workDir ? file(workflow.workDir) : null
    if (workDir?.exists() && !params.debug) {
        workDir.deleteDir()
        }
}

workflow.onError {
    println "Error: ChoCallate execution stopped with the following message: ${workflow.errorMessage}"
}

// Process to create FASTA index file (for FreeBayes)
process CREATE_FAI_INDEX {
    maxForks 1
    cpus 1

    input:
    path(ref_genome)

    output:
    tuple path("${ref_genome}"), path("${ref_genome}.fai"), emit: fai_index

    script:
    """
    samtools faidx --threads ${task.cpus} ${ref_genome}
    """
}

// Process to create sequence dictionary (for GATK4)
process CREATE_SEQ_DICT {
    maxForks 1

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
process BOWTIE2_MAPPING {
    maxForks params.bowtie2_forks
    cpus params.bowtie2_cpu

    tag "${sample_id}"

    input:
    tuple val(sample_id), val(read1), val(read2), val(read3)
    val(ref_index)

    output:
    path("${sample_id}.tmp_bam"), emit: bam

    script:
    if ( params.reads_type == 'pe' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -1 ${read1} -2 ${read2} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.tmp_bam
        """
    else if ( params.reads_type == 'se' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -U ${read1} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.min_map_qual} -F 4 - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.tmp_bam
        """
    else if ( params.reads_type == 'mx' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -1 ${read1} -2 ${read2} -U ${read3} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.bam
        """
    else
        error 'Invalid reads type: ${params.reads_type}. Available types: se, pe, mx'
}

// Process to left-align indels in BAM files
process LEFT_ALIGN_INDELS {
    cpus 1
    maxForks 1

    tag "${bam.baseName}"

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
process INDEXING_BAM {
    maxForks 1
    cpus 1

    tag "${bam.baseName}"

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
process COVERAGE_GENERATION {
    cpus 1
    maxForks 1

    tag "${bam.baseName}"

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
process FREEBAYES {
    maxForks params.freebayes_forks

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    val(ploidy)

    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_freebayes"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_freebayes"), emit: indels_vcf

    script:
    if ( params.reads_source == 'gbs' )
        """
        freebayes --fasta-reference ${ref_genome} --targets ${coverage} --dont-left-align-indels --ploidy ${ploidy} \
            --use-best-n-alleles 4 --min-alternate-qsum ${params.min_base_quality} --hwe-priors-off --no-population-priors \
            --binomial-obs-priors-off --allele-balance-priors-off --min-base-quality ${params.min_base_quality} \
            --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam ${bam} --limit-coverage 250 | \
            bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
            bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
            bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.freebayes

        bcftools view -v snps -Oz -o ${bam.baseName}.snps_freebayes ${bam.baseName}.freebayes

        bcftools view -v indels -Oz -o ${bam.baseName}.indels_freebayes ${bam.baseName}.freebayes
        """

    else if ( params.reads_source == 'wgs' )
        """
        freebayes --fasta-reference ${ref_genome} --targets ${coverage} --dont-left-align-indels --ploidy ${ploidy} \
            --use-best-n-alleles 4 --min-alternate-qsum ${params.min_base_quality} --hwe-priors-off --no-population-priors \
            --allele-balance-priors-off --min-base-quality ${params.min_base_quality} \
            --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam ${bam} --limit-coverage 250 | \
            bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
            bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
            bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.freebayes

        bcftools view -v snps -Oz -o ${bam.baseName}.snps_freebayes ${bam.baseName}.freebayes

        bcftools view -v indels -Oz -o ${bam.baseName}.indels_freebayes ${bam.baseName}.freebayes
        """

    else
        error 'Invalid reads source: ${params.reads_source}. Available sources: gbs, wgs'
}

// Process to call variants using bcftools
process BCFTOOLS {
    maxForks params.bcftools_forks
    cpus params.bcftools_cpu

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    val(ploidy)

    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_bcftools"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_bcftools"), emit: indels_vcf

    script:
    """
    bcftools mpileup -Ou --count-orphans --fasta-ref ${ref_genome} --threads ${task.cpus} --max-depth 250 \
        --min-BQ ${params.min_base_quality} --regions-file ${coverage} ${bam} | \
        bcftools call -Ov --multiallelic-caller --threads ${task.cpus} | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.bcftools

    bcftools view -V indels,mnps,bnd,other -Oz -o ${bam.baseName}.snps_bcftools ${bam.baseName}.bcftools

    bcftools view -v indels -Oz -o ${bam.baseName}.indels_bcftools ${bam.baseName}.bcftools
    """
}

// Process to call variants using GATK4
process GATK4 {
    maxForks params.gatk4_forks

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    path(ref_genome_dict)
    val(ploidy)

    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_gatk"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_gatk"), emit: indels_vcf

    script:
    if (ploidy == 2) {
        """
        gatk HaplotypeCaller -I ${bam} -R ${ref_genome} -mbq ${params.min_base_quality} -O ${bam.baseName}.gatk1 -L ${coverage} --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE
        bcftools filter ${bam.baseName}.gatk1 -e'QUAL<${params.min_snp_qual}' | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools sort - | \
        bcftools view -AA --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.gatk

        bcftools view -v snps -Oz -o ${bam.baseName}.snps_gatk ${bam.baseName}.gatk

        bcftools view -v indels -Oz -o ${bam.baseName}.indels_gatk ${bam.baseName}.gatk
        """
    } else {
        """
        gatk HaplotypeCaller -I ${bam} -R ${ref_genome} -mbq ${params.min_base_quality} -O ${bam.baseName}.gatk1 -L ${coverage} -ERC BP_RESOLUTION -ploidy ${ploidy} --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE
        bcftools filter ${bam.baseName}.gatk1 -e'QUAL<${params.min_snp_qual}' | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools sort - | \
        bcftools view -AA - | bcftools view --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.gatk

        bcftools view -V indels,mnps,bnd,other -Oz -o ${bam.baseName}.snps_gatk ${bam.baseName}.gatk

        bcftools view -v indels -Oz -o ${bam.baseName}.indels_gatk ${bam.baseName}.gatk
        """
    }

}

// Process to call variants using VarDict
process VARDICT {
    maxForks params.vardict_forks
    cpus params.vardict_cpu

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    val(ploidy)

    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_vardict"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_vardict"), emit: indels_vcf

    script:
    """
    vardict-java -G ${ref_genome} -N ${bam.baseName} -b ${bam} -fisher -th ${task.cpus} \
        -VS SILENT --nosv -k 0 -q ${params.min_base_quality} -c 1 -S 2 -E 3 -g 4 ${coverage} | \
        var2vcf_valid.pl -q ${params.min_base_quality} -N ${bam.baseName} -E | \
        bcftools reheader -f ${ref_genome_fai} - | bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.vardict

    bcftools view -v snps -Oz -o ${bam.baseName}.snps_vardict ${bam.baseName}.vardict

    bcftools view -v indels -Oz -o ${bam.baseName}.indels_vardict ${bam.baseName}.vardict
    """
}

// Process to call variants using SNVer
process SNVER {
    maxForks params.snver_forks

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    val(ploidy)

    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_snver"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_snver"), emit: indels_vcf

    script:
    """
    ln -s ${ref_genome} reference.fasta

    samtools faidx reference.fasta

    snver -i ${bam} -r reference.fasta -o ${bam.baseName} -l ${coverage} -bq ${params.min_base_quality} -n ${ploidy}

    bcftools reheader -f reference.fasta.fai ${bam.baseName}.filter.vcf | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Oz -o ${bam.baseName}.snps_snver

    bcftools reheader -f reference.fasta.fai ${bam.baseName}.indel.filter.vcf | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Oz -o ${bam.baseName}.indels_snver
    """
}


// Process to generate final consensus VCF from all callers using majority rule
// Majority rule - variant is true if detected more than 3 callers for diploid calling and 2 for polyploid calling
process GENERATE_CONSENSUS_DIPLOID {
    maxForks params.cons_forks
    cpus params.cons_cpus

    tag "${sample}"

    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.snps.vcf.gz'
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.indels.vcf.gz'

    input:
    tuple val(sample), path(snps_vcf1), path(snps_vcf2), path(snps_vcf3), path(snps_vcf4), path(snps_vcf5)
    tuple val(sample), path(indels_vcf1), path(indels_vcf2), path(indels_vcf3), path(indels_vcf4), path(indels_vcf5)
    path(ref_genome_fai)

    output:
    path("${sample}.snps.vcf.gz")
    path("${sample}.indels.vcf.gz")

    script:
    """
    awk -v OFS='\t' '{print \$1,"0",\$2}' ${ref_genome_fai} > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    for i in ${sample}.snps_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'parallel_cons_diploid.sh {1} {#} ${sample} "snps"' :::: genome_intervals.bed

    bcftools concat --naive-force -Oz all_chrs/*.gz | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.snps.vcf.gz

    rm -r all_chrs/*


    for i in ${sample}.indels_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'parallel_cons_diploid.sh {1} {#} ${sample} "indels"' :::: genome_intervals.bed

    bcftools concat --naive-force -Oz all_chrs/*.gz | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.indels.vcf.gz
    """
}


process GENERATE_CONSENSUS_POLYPLOID {
    maxForks params.cons_forks
    cpus params.cons_cpus

    tag "${sample}"

    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.snps.vcf.gz'
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.indels.vcf.gz'

    input:
    tuple val(sample), path(snps_vcf1), path(snps_vcf2), path(snps_vcf3)
    tuple val(sample), path(indels_vcf1), path(indels_vcf2), path(indels_vcf3)
    path(ref_genome_fai)

    output:
    path("${sample}.snps.vcf.gz")
    path("${sample}.indels.vcf.gz")

    script:
    """
    awk -v OFS='\t' '{print \$1,"0",\$2}' ${ref_genome_fai} > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    for i in ${sample}.snps_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'parallel_cons_polyploid.sh {1} {#} ${sample} "snps"' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.snps.vcf.gz

    rm -r all_chrs/*


    for i in ${sample}.indels_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'parallel_cons_polyploid.sh {1} {#} ${sample} "indels"' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.indels.vcf.gz
    """
}
