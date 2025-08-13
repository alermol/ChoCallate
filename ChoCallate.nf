#!/usr/bin/env nextflow

def available_callers  = 'bcftools,gatk,freebayes,snver,vardict'
def diploid_callers    = 'bcftools,gatk,freebayes,snver,vardict'
def polyploid_callers  = 'gatk,freebayes,snver'
def min_callers_count  = 3

include { getAvailableCallersCount             } from './functions/utils.nf'
include { getConsensusThreshold                } from './functions/utils.nf'
include { allEffectiveCallersInAvailable       } from './functions/utils.nf'
include { effectiveCallersAtLeastThree         } from './functions/utils.nf'
include { allEffectiveCallersDiploidSuitable   } from './functions/utils.nf'
include { allEffectiveCallersPolyploidSuitable } from './functions/utils.nf'

if (!allEffectiveCallersInAvailable(params.effective_callers, available_callers)) {
    exit 1
} else {
    println "All effective callers are in available callers: ${params.effective_callers}"
}

if (effectiveCallersAtLeastThree(params.effective_callers) < min_callers_count) {
    println "The number of effective callers must be at least 3, however, you provided ${getAvailableCallersCount(params.effective_callers)}: ${params.effective_callers}"
    exit 1
}

if (params.ploidy == 2) {
    if (allEffectiveCallersDiploidSuitable(params.effective_callers, diploid_callers)) {
        println "All effective callers are suitable for diploid calling: ${params.effective_callers}"
    } else {
        exit 1
    }
} else if (params.ploidy > 2) {
    if (allEffectiveCallersPolyploidSuitable(params.effective_callers, polyploid_callers)) {
        println "All effective callers are suitable for polyploid calling: ${params.effective_callers}"
    } else {
        exit 1
    }
}

if (params.debug) {
    println "Debug mode is enabled"
} else {
    println "Debug mode is disabled"
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

    GENERATE_ZERO_VCF(INDEXING_BAM.out.ind_bam,
                      CREATE_FAI_INDEX.out.fai_index,
                      COVERAGE_GENERATION.out.coverage)

    CALLING(INDEXING_BAM.out.ind_bam,
            CREATE_FAI_INDEX.out.fai_index,
            COVERAGE_GENERATION.out.coverage,
            params.ploidy,
            CREATE_SEQ_DICT.out.gen_dict)

    // Generate final consensus
    if (params.ploidy == 2) {
        GENERATE_CONSENSUS_DIPLOID(CALLING.out.merged_snps_vcfs,
                                   CALLING.out.merged_indels_vcfs,
                                   CREATE_FAI_INDEX.out.fai_index.map{it[1]},
                                   GENERATE_ZERO_VCF.out.zero_vcf)
    } else {
        GENERATE_CONSENSUS_POLYPLOID(CALLING.out.merged_snps_vcfs,
                                     CALLING.out.merged_indels_vcfs,
                                     CREATE_FAI_INDEX.out.fai_index.map{it[1]},
                                     GENERATE_ZERO_VCF.out.zero_vcf)
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


process GENERATE_ZERO_VCF {
    maxForks params.zero_vcf_forks
    cpus params.zero_vcf_cpu

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    tuple path(ref_genome), path(ref_genome_fai)
    path(coverage_bed)

    output:
    path("zero.vcf.gz"), emit: zero_vcf

    script:
    String genotype = (["0"] * params.ploidy).join("/")
    """
    bcftools mpileup -Ov --count-orphans --fasta-ref ${ref_genome} --threads ${task.cpus} --max-depth 1 \
        --min-BQ ${params.min_base_quality} --regions-file ${coverage_bed} ${bam} | \
        awk -v OFS='\\t' -v gen=${genotype} '{if(\$0 !~ /#/) print \$1,\$2,\$3,\$4,".","100",".",".","GT",gen; else print \$0}' | \
        awk -v OFS='\\t' '{if(length(\$4) == 1 || \$0 ~ /#/) print \$0}' | bgzip  > zero.vcf.gz
    """
}

process CALLING {
    maxForks 1
    cpus 1

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    tuple path(ref_genome), path(ref_genome_fai)
    path(coverage_bed)
    val(ploidy)
    path(ref_genome_dict)

    output:
    tuple val("${bam.baseName}"), path("*.snps_*.vcf.gz"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("*.indels_*.vcf.gz"), emit: indels_vcf

    script:
    def parallel_cpus = getAvailableCallersCount(params.effective_callers)
    """
    touch callers_commands.sh

    if [[ ",${params.effective_callers}," == *"bcftools"* ]]; then
        echo "bash ${projectDir}/bin/bcftools_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.bcftools_cpu}" >> callers_commands.sh
    fi

    if [[ ",${params.effective_callers}," == *"freebayes"* ]]; then
        echo "bash ${projectDir}/bin/freebayes_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.reads_source} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
    fi

    if [[ ",${params.effective_callers}," == *"gatk"* ]]; then
        echo "bash ${projectDir}/bin/gatk4_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
    fi

    if [[ ",${params.effective_callers}," == *"vardict"* ]]; then
        echo "bash ${projectDir}/bin/vardict_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.vardict_cpu}" >> callers_commands.sh
    fi

    if [[ ",${params.effective_callers}," == *"snver"* ]]; then
        echo "bash ${projectDir}/bin/snver_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
    fi

    parallel -j ${parallel_cpus} '{}' :::: callers_commands.sh
    """
}


// Process to generate final consensus VCF
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
    path(zero_vcf)

    output:
    path("${sample}.snps.vcf.gz")
    path("${sample}.indels.vcf.gz")

    script:
    def cons_threshold = getConsensusThreshold(params.cons_type, available_callers)
    """
    awk -v OFS='\t' '{print \$1,"0",\$2}' ${ref_genome_fai} > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    for i in ${sample}.snps_*; do tabix -C \${i}; done

    tabix -C zero.vcf.gz

    parallel -j ${task.cpus} 'parallel_cons_diploid.sh {1} {#} ${sample} "snps" ${cons_threshold}' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.snps.vcf.gz

    rm -r all_chrs/*


    for i in ${sample}.indels_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'parallel_cons_diploid.sh {1} {#} ${sample} "indels" ${cons_threshold}' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
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
    path(zero_vcf)

    output:
    path("${sample}.snps.vcf.gz")
    path("${sample}.indels.vcf.gz")

    script:
    def cons_threshold = getConsensusThreshold(params.cons_type, available_callers)
    """
    awk -v OFS='\t' '{print \$1,"0",\$2}' ${ref_genome_fai} > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    for i in ${sample}.snps_*; do tabix -C \${i}; done

    tabix -C zero.vcf.gz

    parallel -j ${task.cpus} 'parallel_cons_polyploid.sh {1} {#} ${sample} "snps" ${cons_threshold}' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.snps.vcf.gz

    rm -r all_chrs/*


    for i in ${sample}.indels_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'parallel_cons_polyploid.sh {1} {#} ${sample} "indels" ${cons_threshold}' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.indels.vcf.gz
    """
}
