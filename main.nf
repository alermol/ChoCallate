#!/usr/bin/env nextflow

// Generate sample channel
include { GENERATE_SAMPLE_CHANNEL } from './modules/generate_channel.nf'

// Prepare assembly
include { DECOMPRESS_ASSEMBLY } from './modules/decompress_assembly.nf'
include { CREATE_FAI_INDEX } from './modules/create_fai_index.nf'
include { CREATE_SEQ_DICT } from './modules/create_seq_dict.nf'

// Map reads to assembly using Bowtie2
include { MAP_BOWTIE2_SINGLE; MAP_BOWTIE2_PAIRED; MAP_BOWTIE2_MIXED } from './modules/mapping_bowtie2.nf'

// Map reads to assembly using BWA
include { MAP_BWA_SINGLE; MAP_BWA_PAIRED } from './modules/mapping_bwa.nf'

// Map reads to assembly using Minimap2
include { MAP_MINIMAP2_SINGLE; MAP_MINIMAP2_PAIRED } from './modules/mapping_minimap2.nf'

// Remove duplicates
include { REMOVE_DUPLICATES } from './modules/remove_duplicates.nf'

// Filter BAM file
include { FILTER_MAPPING_BAM; FILTER_INPUT_BAM } from './modules/filter_bam.nf'

// Left align indels
include { LEFT_ALIGN_INDELS } from './modules/left_align_indels.nf'

// Generate coverage
include { GENERATE_COVERAGE; INTERSECT_CUSTOM_BED } from './modules/generate_coverage.nf'

// Call variants
include { CALL_BCFTOOLS; CALL_FREEBAYES; CALL_GATK; CALL_SNVER; CALL_VARDICT } from './modules/calling.nf'

// Generate consensus
include { GENERATE_CONSENSUS } from './modules/generate_consensus.nf'

// Convert consensus to VCF
include { CONVERT_TO_VCF } from './modules/convert_to_vcf.nf'

// Merge outputs
include { MERGE_OUTPUTS } from './modules/merge_outputs.nf'


workflow {

    // Validate required parameters
    CLIParamsValidation.outdir_validation(params.outdir)
    CLIParamsValidation.cons_threshold_validation(params.cons_threshold, params.callers)
    CLIParamsValidation.reference_genome_validation(params.reference_genome)
    CLIParamsValidation.reference_index_validation(params.reference_index, params.input_format)
    CLIParamsValidation.samples_tsv_validation(params.samples_tsv)
    CLIParamsValidation.ploidy_validation(params.ploidy, workflow.profile)
    CLIParamsValidation.effective_callers_validation(params.callers, workflow.profile, params.ploidy, params.cons_threshold)
    CLIParamsValidation.mapper_validation(params.bam_preparation.mapping.mapper, params.reads_type)

    sample_run_ch = GENERATE_SAMPLE_CHANNEL(params.samples_tsv, params.input_format, params.reads_type)

    ref_genome = file(params.reference_genome)

    if (params.reference_preparation.decompress.bgzip) {
        DECOMPRESS_ASSEMBLY(ref_genome)
        ref_genome = DECOMPRESS_ASSEMBLY.out.ref_genome
    }

    fai_index = CREATE_FAI_INDEX(ref_genome)
    gen_dict = CREATE_SEQ_DICT(ref_genome)

    if (params.input_format == 'fastq') {
        ref_index = file(params.reference_index)
        if (params.bam_preparation.mapping.mapper == 'bowtie2') {
            if (params.reads_type == 'pe') {
                bam_file = MAP_BOWTIE2_PAIRED(sample_run_ch, ref_index)
            } else if (params.reads_type == 'se') {
                bam_file = MAP_BOWTIE2_SINGLE(sample_run_ch, ref_index)
            } else if (params.reads_type == 'mx') {
                bam_file = MAP_BOWTIE2_MIXED(sample_run_ch, ref_index)
            }
        }
        if (params.bam_preparation.mapping.mapper == 'bwa') {
            if (params.reads_type == 'pe') {
                bam_file = MAP_BWA_PAIRED(sample_run_ch, ref_index)
            } else if (params.reads_type == 'se') {
                bam_file = MAP_BWA_SINGLE(sample_run_ch, ref_index)
            }
        }
        if (params.bam_preparation.mapping.mapper == 'minimap2') {
            if (params.reads_type == 'pe') {
                bam_file = MAP_MINIMAP2_PAIRED(sample_run_ch, ref_index)
            } else if (params.reads_type == 'se') {
                bam_file = MAP_MINIMAP2_SINGLE(sample_run_ch, ref_index)
            }
        }
        bam_file = FILTER_MAPPING_BAM(bam_file)
    }

    if (params.input_format == 'bam') {
        bam_file = FILTER_INPUT_BAM(sample_run_ch)
    }

    if (params.bam_preparation.rm_duplicates.enabled) {
        bam_file = REMOVE_DUPLICATES(bam_file)
    }

    if (params.bam_preparation.left_align_indels.enabled) {
        bam_file = LEFT_ALIGN_INDELS(bam_file, ref_genome, gen_dict, fai_index)
    }

    bed_coverage = GENERATE_COVERAGE(bam_file)
    if (params.custom_bed != null) {
        bed_coverage = INTERSECT_CUSTOM_BED(bed_coverage, file(params.custom_bed))
    }

    if (params.callers.contains('bcftools')) {
        CALL_BCFTOOLS(bam_file, ref_genome, bed_coverage)
        bcftools = CALL_BCFTOOLS.out.calling_result
    } else {
        bcftools = channel.empty()
    }

    if (params.callers.contains('freebayes')) {
        CALL_FREEBAYES(bam_file, ref_genome, bed_coverage)
        freebayes = CALL_FREEBAYES.out.calling_result
    } else {
        freebayes = channel.empty()
    }

    if (params.callers.contains('gatk')) {
        CALL_GATK(bam_file, ref_genome, fai_index, bed_coverage, gen_dict)
        gatk = CALL_GATK.out.calling_result
    } else {
        gatk = channel.empty()
    }

    if (params.callers.contains('snver')) {
        CALL_SNVER(bam_file, ref_genome, fai_index, bed_coverage)
        snver = CALL_SNVER.out.calling_result
    } else {
        snver = channel.empty()
    }

    if (params.callers.contains('vardict')) {
        CALL_VARDICT(bam_file, ref_genome, fai_index, bed_coverage)
        vardict = CALL_VARDICT.out.calling_result
    } else {
        vardict = channel.empty()
    }

    all_calls = bcftools
        .join(freebayes, remainder: true)
        .join(gatk, remainder: true)
        .join(snver, remainder: true)
        .join(vardict, remainder: true)
        .map { list -> list.findAll { item -> item != null } }
        .map {tuple -> [tuple[0], tuple[1..-1]]}



    // Generate BCF consensus file for each sample
    if (params.output.format == 'bcf' && params.output.type == 'sample') {
        GENERATE_CONSENSUS(all_calls, ref_genome, fai_index, bed_coverage)
    }
    
    // Generate single BCF consensus file for all samples
    if (params.output.format == 'bcf' && params.output.type == 'single') {
        GENERATE_CONSENSUS(all_calls, ref_genome, fai_index, bed_coverage)
        consensus = GENERATE_CONSENSUS.out.consensus.map { item -> item[1] }.collect()
        MERGE_OUTPUTS(consensus)
    }
    
    // Generate VCF consensus file for each sample
    if (params.output.format == 'vcf' && params.output.type == 'sample') {
        GENERATE_CONSENSUS(all_calls, ref_genome, fai_index, bed_coverage)
        CONVERT_TO_VCF(GENERATE_CONSENSUS.out.consensus)
    }
    
    // Generate single VCF consensus file for all samples
    if (params.output.format == 'vcf' && params.output.type == 'single') {
        GENERATE_CONSENSUS(all_calls, ref_genome, fai_index, bed_coverage)
        consensus = GENERATE_CONSENSUS.out.consensus.map { item -> item[1] }.collect()
        consensus = MERGE_OUTPUTS(consensus)
        CONVERT_TO_VCF(consensus)
    }
    
}

