#!/usr/bin/env nextflow

include { DECOMPRESS_ASSEMBLY } from './modules/local/decompress_assembly'
include { CREATE_FAI_INDEX } from './modules/local/create_fai_index'
include { CREATE_SEQ_DICT } from './modules/local/create_seq_dict'
include { MAPPING_BOWTIE2 } from './modules/local/mapping_bowtie2'
include { MAPPING_BWA } from './modules/local/mapping_bwa'
include { MAPPING_MINIMAP2 } from './modules/local/mapping_minimap2'
include { PREPARE_BAM } from './modules/local/prepare_bam'
include { GENERATE_COVERAGE } from './modules/local/generate_coverage'
include { SCATTER_INTERVALS } from './modules/local/scatter_intervals'
include { CALLING_BCFTOOLS } from './modules/local/calling_bcftools'
include { CALLING_FREEBAYES } from './modules/local/calling_freebayes'
include { CALLING_GATK } from './modules/local/calling_gatk'
include { GENERATE_CONSENSUS } from './modules/local/generate_consensus'
include { MERGE_OUTPUTS } from './modules/local/merge_outputs'
include { validateParameters; paramsHelp } from 'plugin/nf-schema'

workflow {

    if (params.help || params.helpFull) {
        def help_options = [command: "nextflow run main.nf --help"]
        def help_value = params.help?.toString() ?: 'false'
        if (help_value != 'true' && help_value != 'false') {
            help_options['parameter'] = help_value
        }
        if (params.helpFull?.toString() == 'true') {
            help_options['fullHelp'] = true
        }
        // nf-schema filters its built-in help options by line prefix; the
        // wrapped continuation lines of those descriptions would otherwise
        // remain in the output.
        def help_text = paramsHelp(help_options)
            .readLines()
            .findAll { line ->
                !line.startsWith(
                    'help message of that parameter will be printed.'
                ) && !line.startsWith('or `--helpFull`.')
            }
            .join('\n')
        log.info help_text
        exit 0
    }

    // Unknown keys are rejected before any process starts.
    validateParameters()

    // Caller selection: trimmed, deduplicated, in voting order. The same list
    // drives process selection, threshold validation and the consensus header.
    def selected_callers =
        CLIParamsValidation.parse_callers(params.calling.callers)
    CLIParamsValidation.effective_callers_validation(params.calling.callers)
    CLIParamsValidation.cons_threshold_validation(
        params.consensus.threshold,
        params.calling.callers
    )
    CLIParamsValidation.mixed_reads_validation(
        params.mapping.mapper,
        params.input.reads_type,
        params.input.format
    )

    // Parse the samplesheet up front: unique shell-safe ids, exact row widths,
    // optional header/comment lines, and read paths relative to the sheet.
    def samples = CLIParamsValidation.parse_samplesheet(
        file(params.input.samples_tsv, checkIfExists: true),
        params.input.format,
        params.input.reads_type
    )

    // A merge opens every sample BCF plus the named pipes; fail before the
    // first task when the process descriptor limit cannot support that.
    if (params.output.type == 'single') {
        CLIParamsValidation.merge_ulimit_validation(
            samples.size(),
            params.consensus.cpu as int
        )
    }

    // Create the sample channel.  Absent read slots stage the empty NO_FILE
    // instead of /dev/null so every input is a real, portable file.
    def no_file = file("${projectDir}/assets/NO_FILE", checkIfExists: true)
    sample_run_ch = channel.fromList(
        samples.collect { sample ->
            def reads = sample.reads.collect { read_path -> file(read_path) }
            reads = reads + ([no_file] * (3 - reads.size()))
            tuple(sample.id, reads[0], reads[1], reads[2])
        }
    )

    // Prepare genome assembly, fai index and sequence dictionary.  BGZF is
    // detected from the header, so ref_genome.bgzip is no longer a parameter.
    def reference_file =
        file(params.input.reference_genome, checkIfExists: true)
    def reference_is_bgzip = CLIParamsValidation.is_bgzip(reference_file)
    ref_genome = reference_is_bgzip
        ? DECOMPRESS_ASSEMBLY(reference_file)
        : reference_file
    fai_index = CREATE_FAI_INDEX(ref_genome)
    gen_dict = CREATE_SEQ_DICT(ref_genome)

    // Prepare BAM files
    if (params.mapping.mapper == 'bowtie2' && params.input.format != 'bam') {
        bam_file = MAPPING_BOWTIE2(sample_run_ch,
                              file(params.input.reference_genome),
                              ref_genome,
                              gen_dict,
                              fai_index)
    } else if (params.mapping.mapper == 'bwa' &&
               params.input.format != 'bam') {
        bam_file = MAPPING_BWA(sample_run_ch, 
                          file(params.input.reference_genome), 
                          ref_genome, 
                          gen_dict, 
                          fai_index)  
    } else if (params.mapping.mapper == 'minimap2' &&
               params.input.format != 'bam') {
        bam_file = MAPPING_MINIMAP2(sample_run_ch, 
                               file(params.input.reference_genome), 
                               ref_genome, 
                               gen_dict, 
                               fai_index)
    } else if (params.input.format == 'bam') {
        bam_file = PREPARE_BAM(sample_run_ch,
                          ref_genome, 
                          gen_dict, 
                          fai_index)
    } else {
        error("Unknown input format or mapper")
    }

    // Generate BED coverage
    include_bed = params.input.include_bed == null
        ? file("${projectDir}/assets/N1_FILE", checkIfExists: true)
        : file(params.input.include_bed, checkIfExists: true)
    exclude_bed = params.input.exclude_bed == null
        ? file("${projectDir}/assets/NO_FILE", checkIfExists: true)
        : file(params.input.exclude_bed, checkIfExists: true)
    bed_coverage = GENERATE_COVERAGE(bam_file, include_bed, exclude_bed)

    // Scatter the covered intervals once per sample. Every caller and the
    // consensus step consume the same chunk list instead of rebuilding it.
    intervals = SCATTER_INTERVALS(bed_coverage, ref_genome, fai_index, gen_dict)

    // Join each sample's BAM, BAM index and interval chunks by sample_id. This
    // keeps all per-sample inputs coupled and avoids positional mixing.
    bam_intervals = bam_file.join(intervals)

    // Perform variant calling
    bcftools = selected_callers.contains('bcftools')
        ? CALLING_BCFTOOLS(bam_intervals, ref_genome, fai_index, gen_dict)
        : channel.empty()
    freebayes = selected_callers.contains('freebayes')
        ? CALLING_FREEBAYES(bam_intervals, ref_genome, fai_index, gen_dict)
        : channel.empty()
    gatk = selected_callers.contains('gatk')
        ? CALLING_GATK(bam_intervals, ref_genome, fai_index, gen_dict)
        : channel.empty()
    all_calls = bcftools
        .join(freebayes, remainder: true)
        .join(gatk, remainder: true)
        .map { list -> list.findAll { item -> item != null } }
        .map { tuple ->
            def paths = tuple[1..-1]
            def received = paths.collect { path ->
                path.name.replaceAll(/\.bcf$/, '')
            }
            def failed = selected_callers.findAll { caller ->
                !(caller in received)
            }
            if (failed) {
                log.warn(
                    "Sample ${tuple[0]}: caller(s) failed and were " +
                    "skipped: ${failed.join(',')}"
                )
            }
            [tuple[0], received.join(','), paths]
        }

    // Generate consensus.  The configured caller list is passed separately so
    // the module can record the structured consensus header and warn about
    // degrades.
    GENERATE_CONSENSUS(
        all_calls.join(bam_intervals),
        selected_callers.join(','),
        ref_genome,
        fai_index
    )

    // Merge consensuses from different samples into single VCF of BCF
    if (params.output.type == 'single') {
        consensus_ch = params.output.format == 'vcf' ?
                       GENERATE_CONSENSUS.out.consensus_vcf :
                       GENERATE_CONSENSUS.out.consensus_bcf
        sample_rank = sample_run_ch
                        .map { item -> item[0] }
                        .collect()
                        .flatMap { ids ->
                            ids.withIndex().collect { id, i ->
                                tuple(id, i)
                            }
                        }
        consensus = consensus_ch
                        .join(sample_rank)
                        .toSortedList { item -> item[2] }
                        .map { items ->
                            def ids = items.collect { item -> item[0] }
                            [ids.join(','), items.collect { item -> item[1] }]
                        }
        MERGE_OUTPUTS(consensus, ref_genome, fai_index)
    }
}

