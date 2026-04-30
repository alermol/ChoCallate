#!/usr/bin/env nextflow

include { DECOMPRESS_ASSEMBLY } from './modules/local/decompress_assembly'
include { CREATE_FAI_INDEX } from './modules/local/create_fai_index'
include { CREATE_SEQ_DICT } from './modules/local/create_seq_dict'
include { MAPPING_BOWTIE2 } from './modules/local/mapping_bowtie2'
include { MAPPING_BWA } from './modules/local/mapping_bwa'
include { MAPPING_MINIMAP2 } from './modules/local/mapping_minimap2'
include { PREPARE_BAM } from './modules/local/prepare_bam'
include { GENERATE_COVERAGE } from './modules/local/generate_coverage'
include { CALLING_BCFTOOLS } from './modules/local/calling_bcftools'
include { CALLING_FREEBAYES } from './modules/local/calling_freebayes'
include { CALLING_GATK } from './modules/local/calling_gatk'
include { GENERATE_CONSENSUS } from './modules/local/generate_consensus'
include { MERGE_OUTPUTS } from './modules/local/merge_outputs'
include { MERGE_BEDS } from './modules/local/merge_beds'


workflow {
    
    // Validate required parameters
    CLIParamsValidation.reference_genome_validation(params.input.reference_genome)
    CLIParamsValidation.samples_tsv_validation(params.input.samples_tsv)
    CLIParamsValidation.effective_callers_validation(params.calling.callers)
    CLIParamsValidation.cons_threshold_validation(params.consensus.threshold, params.calling.callers)
    CLIParamsValidation.mapper_validation(params.mapping.mapper, params.input.reads_type)
    CLIParamsValidation.reference_index_dir_validation(params.input.reference_index_dir)

    // Create sample channel
    sample_run_ch = channel
                .fromPath(params.input.samples_tsv)
                .splitCsv(header: false, sep: '\t')
                .filter { row ->
                    if (row.size() == 1) {
                        println "[WARN] Dropping row with only one column: ${row[0]}"
                        return false
                    }
                    return true
                }
                .map { 
                    row -> if (params.input.input_format == 'bam' || params.input.reads_type == 'se') {
                        tuple(row[0..1].toArray() + ['/dev/null'] * 2)
                    } else if (params.input.reads_type == 'pe' && params.input.format != 'bam') {
                        tuple(row[0..2].toArray() + ['/dev/null'])
                    } else if (params.input.reads_type == 'mx' && params.input.format != 'bam') {
                        tuple(row[0..3].toArray())
                    }
                 }

    // Prepare genome assembly, fai index and sequence dictionary
    ref_genome = params.ref_genome.bgzip ? DECOMPRESS_ASSEMBLY(file(params.input.reference_genome)) : file(params.input.reference_genome)
    fai_index = CREATE_FAI_INDEX(ref_genome)
    gen_dict = CREATE_SEQ_DICT(ref_genome)

    // Prepare BAM files
    bam_file = params.mapping.mapper == 'bowtie2' && params.input.format != 'bam' ? MAPPING_BOWTIE2(sample_run_ch, 
                                                                                                    file(params.input.reference_genome), 
                                                                                                    ref_genome, 
                                                                                                    gen_dict, 
                                                                                                    fai_index) : 
               params.mapping.mapper == 'bwa' && params.input.format != 'bam' ? MAPPING_BWA(sample_run_ch, 
                                                                                            file(params.input.reference_genome), 
                                                                                            ref_genome, 
                                                                                            gen_dict, 
                                                                                            fai_index) : 
               params.mapping.mapper == 'minimap2' && params.input.format != 'bam' ? MAPPING_MINIMAP2(sample_run_ch, 
                                                                                                      file(params.input.reference_genome), 
                                                                                                      ref_genome, 
                                                                                                      gen_dict, 
                                                                                                      fai_index) : 
               params.input.format == 'bam' ? PREPARE_BAM(sample_run_ch, 
                                                          file(params.input.reference_genome), 
                                                          ref_genome, 
                                                          gen_dict, 
                                                          fai_index) : 
               error("Unknown input format or mapper")

    // Generate BED coverage
    include_bed = params.input.include_bed == null ? file("${projectDir}/assets/N1_FILE", checkIfExists: true) : file(params.input.include_bed, checkIfExists: true)
    exclude_bed = params.input.exclude_bed == null ? file("${projectDir}/assets/NO_FILE", checkIfExists: true) : file(params.input.exclude_bed, checkIfExists: true)
    bed_coverage = GENERATE_COVERAGE(bam_file, include_bed, exclude_bed)

    // Perform variant calling
    bcftools = params.calling.callers.contains('bcftools') ? CALLING_BCFTOOLS(bam_file, ref_genome, fai_index, bed_coverage) : channel.empty()
    freebayes = params.calling.callers.contains('freebayes') ? CALLING_FREEBAYES(bam_file, ref_genome, fai_index, bed_coverage) : channel.empty()
    gatk = params.calling.callers.contains('gatk') ? CALLING_GATK(bam_file, ref_genome, fai_index, bed_coverage, gen_dict) : channel.empty()
    all_calls = bcftools
        .join(freebayes, remainder: true)
        .join(gatk, remainder: true)
        .map { list -> list.findAll { item -> item != null } }
        .map {tuple -> [tuple[0], tuple[1..-1]]}

    // Generate consensus
    GENERATE_CONSENSUS(all_calls, ref_genome, fai_index, bed_coverage)

    // Merge consensuses from different samples into single VCF of BCF
    if (params.output.type == 'single') {
        consensus = params.output.format == 'vcf' ? 
                    GENERATE_CONSENSUS.out.consensus_vcf.map { item -> item[1] }.collect() : 
                    GENERATE_CONSENSUS.out.consensus_bcf.map { item -> item[1] }.collect()
        bed_coverage = MERGE_BEDS(bed_coverage.collect())
        MERGE_OUTPUTS(consensus, bed_coverage)
    }
}

