#!/usr/bin/env nextflow

include { DECOMPRESS_ASSEMBLY } from './modules/local/decompress_assembly'
include { CREATE_FAI_INDEX } from './modules/local/create_fai_index'
include { CREATE_SEQ_DICT } from './modules/local/create_seq_dict'
include { MAPPING_BOWTIE2 } from './modules/local/mapping_bowtie2'
include { MAPPING_BWA } from './modules/local/mapping_bwa'
include { MAPPING_MINIMAP2 } from './modules/local/mapping_minimap2'
include { MAPPING_STATS } from './modules/local/mapping_stats'
include { MAPPING_SUMMARY } from './modules/local/mapping_summary'
include { PREPARE_BAM } from './modules/local/prepare_bam'
include { GENERATE_COVERAGE } from './modules/local/generate_coverage'
include { SCATTER_INTERVALS } from './modules/local/scatter_intervals'
include { CALLING_BCFTOOLS } from './modules/local/calling_bcftools'
include { CALLING_FREEBAYES } from './modules/local/calling_freebayes'
include { CALLING_GATK } from './modules/local/calling_gatk'
include { GENERATE_CONSENSUS } from './modules/local/generate_consensus'
include { MERGE_OUTPUTS } from './modules/local/merge_outputs'
include { REPORT_DROPPED_SAMPLES } from './modules/local/report_dropped_samples'
include { EXPORT_OUTPUT } from './modules/local/export_output'
include { validateParameters } from 'plugin/nf-schema'

workflow {

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
    CLIParamsValidation.ploidy_validation(params.calling)
    CLIParamsValidation.read_group_validation(params.mapping)

    // Parse the samplesheet up front: unique shell-safe ids, exact row widths,
    // optional header/comment lines, and read paths relative to the sheet.
    def samples = CLIParamsValidation.parse_samplesheet(
        file(params.input.samples_tsv, checkIfExists: true),
        params.input.format,
        params.input.reads_type
    )
    def sample_ids = samples.collect { sample -> sample.id }

    // A merge opens every sample BCF plus the named pipes; fail before the
    // first task when the process descriptor limit cannot support that.
    if (params.output.type == 'single') {
        CLIParamsValidation.merge_ulimit_validation(
            samples.size(),
            params.consensus.cpu as int
        )
    }

    // Compression is detected from the header.  BGZF is decompressed below;
    // plain gzip cannot be indexed by samtools faidx or read by GATK, so it
    // is rejected before any task starts.
    def reference_file =
        file(params.input.reference_genome, checkIfExists: true)
    def reference_is_bgzip = CLIParamsValidation.is_bgzip(reference_file)
    if (!reference_is_bgzip && CLIParamsValidation.is_gzip(reference_file)) {
        error(
            "reference is gzip-compressed but not BGZF; recompress with bgzip"
        )
    }

    // Mapper indexes are staged as task inputs. The prefix is resolved here,
    // so the command never depends on a host path, and a missing index aborts
    // before the first task starts. Exact suffixes keep the FASTA, its
    // .fai/.dict and indexes of other genomes that share the prefix out of
    // the staged set. An index is usable only when every file of one set is
    // present.
    def index_suffix_sets = [
        bowtie2: [
            ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2',
             '.rev.1.bt2', '.rev.2.bt2'],
            ['.1.bt2l', '.2.bt2l', '.3.bt2l', '.4.bt2l',
             '.rev.1.bt2l', '.rev.2.bt2l']
        ],
        bwa: [
            ['.amb', '.ann', '.bwt', '.pac', '.sa']
        ],
        minimap2: [
            ['.mmi']
        ]
    ]
    def index_files = channel.empty()
    if (params.input.format != 'bam') {
        def mapper = params.mapping.mapper
        def index_dir = params.input.reference_index_dir
        CLIParamsValidation.reference_index_prefix_validation(
            params.input.reference_index_prefix
        )
        def index_prefix = params.input.reference_index_prefix ?:
            file(params.input.reference_genome).name
        def suffix_sets = index_suffix_sets[mapper]
        if (suffix_sets == null) {
            error("Unknown mapper '${mapper}'")
        }
        def candidate_sets = suffix_sets.collect { suffixes ->
            suffixes.collect { suffix ->
                file("${index_dir}/${index_prefix}${suffix}")
            }
        }
        def index_matches = candidate_sets.find { candidates ->
            candidates.every { candidate -> candidate.exists() }
        }
        if (index_matches == null) {
            def missing = candidate_sets.collect { candidates ->
                candidates
                    .findAll { candidate -> !candidate.exists() }
                    .collect { candidate -> candidate.name }
                    .join(', ')
            }
            error(
                "No ${mapper} index files found in '${index_dir}' for " +
                "prefix '${index_prefix}'; missing: ${missing.join('; ')}"
            )
        }
        index_files = channel.value(index_matches)
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

    // Prepare genome assembly, fai index and sequence dictionary.
    ref_genome = reference_is_bgzip
        ? DECOMPRESS_ASSEMBLY(reference_file)
        : reference_file
    fai_index = CREATE_FAI_INDEX(ref_genome)
    gen_dict = CREATE_SEQ_DICT(ref_genome)

    // Prepare BAM files
    if (params.mapping.mapper == 'bowtie2' && params.input.format != 'bam') {
        bam_file = MAPPING_BOWTIE2(sample_run_ch,
                              index_files,
                              ref_genome,
                              gen_dict,
                              fai_index)
    } else if (params.mapping.mapper == 'bwa' &&
               params.input.format != 'bam') {
        bam_file = MAPPING_BWA(sample_run_ch,
                          index_files,
                          ref_genome,
                          gen_dict,
                          fai_index)
    } else if (params.mapping.mapper == 'minimap2' &&
               params.input.format != 'bam') {
        bam_file = MAPPING_MINIMAP2(sample_run_ch,
                               index_files,
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

    // The raw flagstat is captured by the mapping stage before MAPQ filtering
    // and duplicate handling; the other reports describe the final BAM.
    mapping_stats = MAPPING_STATS(bam_file)
    raw_flagstat_files = mapping_stats
        .map { _sample_id, raw, _flagstat, _stats, _idxstats -> raw }
        .collect()
    final_flagstat_files = mapping_stats
        .map { _sample_id, _raw, flagstat, _stats, _idxstats -> flagstat }
        .collect()
    mapping_stats_files = mapping_stats
        .map { _sample_id, _raw, _flagstat, stats, _idxstats -> stats }
        .collect()
    mapping_idxstats_files = mapping_stats
        .map { _sample_id, _raw, _flagstat, _stats, idxstats -> idxstats }
        .collect()
    MAPPING_SUMMARY(
        channel.value(sample_ids),
        raw_flagstat_files,
        final_flagstat_files,
        mapping_stats_files,
        mapping_idxstats_files
    )

    // Coverage, calling and consensus only need the BAM and its index.
    bam_for_analysis = bam_file.map { sample_id, bam, csi, _raw ->
        tuple(sample_id, bam, csi)
    }

    // Generate BED coverage
    include_bed = params.input.include_bed == null
        ? file("${projectDir}/assets/N1_FILE", checkIfExists: true)
        : file(params.input.include_bed, checkIfExists: true)
    exclude_bed = params.input.exclude_bed == null
        ? file("${projectDir}/assets/NO_FILE", checkIfExists: true)
        : file(params.input.exclude_bed, checkIfExists: true)
    bed_coverage = GENERATE_COVERAGE(
        bam_for_analysis, include_bed, exclude_bed, fai_index
    )

    // With coverage.on_empty="skip", an empty coverage BED removes the
    // sample before scattering.  The default "fail" mode is handled inside
    // GENERATE_COVERAGE, which aborts with the sample id and threshold.
    if (params.coverage.on_empty == 'skip') {
        coverage_branch = bed_coverage.branch { _sample_id, bed ->
            empty: bed.size() == 0
            covered: true
        }
        empty_coverage_samples = coverage_branch.empty
            .map { sample_id, _bed ->
                log.warn(
                    "Sample ${sample_id}: no region at or above " +
                    "coverage.min_coverage=" +
                    "${params.coverage.min_coverage}; sample dropped"
                )
                sample_id
            }
            .toList()
        bed_coverage = coverage_branch.covered
    } else {
        empty_coverage_samples = channel.value([])
    }

    // Scatter the covered intervals once per sample. Every caller and the
    // consensus step consume the same chunk list instead of rebuilding it.
    intervals = SCATTER_INTERVALS(
        bed_coverage, ref_genome, fai_index, gen_dict
    )

    // Join each sample's BAM, BAM index and interval chunks by sample_id. This
    // keeps all per-sample inputs coupled and avoids positional mixing.
    bam_intervals = bam_for_analysis.join(intervals)

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

    // A sample whose every caller failed produces no joined tuple and would
    // otherwise disappear from the output without a message.  Compare the ids
    // that reached the join with the ids that passed the coverage filter.
    coverage_ids = bed_coverage.map { item -> item[0] }.toList()
    call_ids = all_calls.map { item -> item[0] }.toList()
    dropped_samples = coverage_ids
        .concat(empty_coverage_samples)
        .concat(call_ids)
        .toList()
        .map { lists ->
            def covered = lists[0]
            def empty = lists[1]
            def called = lists[2]
            def failed = covered.findAll { id -> !(id in called) }
            failed.each { id ->
                def message = "Sample ${id}: no caller produced calls " +
                    "(${selected_callers.join(',')} configured); " +
                    "sample dropped from consensus"
                if (params.calling.allow_caller_failure) {
                    log.warn(message)
                } else {
                    error(message)
                }
            }
            (empty + failed).unique().sort()
        }

    // Always publish the dropped list, even when it is empty. The output is
    // used below to order the all-dropped error after the file is published.
    dropped_report = REPORT_DROPPED_SAMPLES(dropped_samples)

    // Generate consensus.  The per-sample consensus is always an internal BCF.
    // Reference records are written whenever a later merge may need them;
    // per-sample output that discards them sets --variants_only instead.
    GENERATE_CONSENSUS(
        all_calls.join(bam_intervals),
        selected_callers.join(','),
        ref_genome,
        fai_index
    )

    // A single output is one multi-sample file: merge the internal BCFs first.
    // Both paths end in EXPORT_OUTPUT, which applies output.format, output.type
    // and the final record transformations.
    def all_dropped_message =
        "All samples were dropped; see dropped_samples.txt"
    if (params.output.type == 'single') {
        sample_rank = sample_run_ch
                        .map { item -> item[0] }
                        .collect()
                        .flatMap { ids ->
                            ids.withIndex().collect { id, i ->
                                tuple(id, i)
                            }
                        }
        consensus = GENERATE_CONSENSUS.out.consensus_bcf
                        .combine(dropped_report)
                        .map { sample_id, bcf, _report ->
                            tuple(sample_id, bcf)
                        }
                        .join(sample_rank)
                        .toSortedList { item -> item[2] }
                        .map { items ->
                            def ids = items.collect { item -> item[0] }
                            [ids.join(','), items.collect { item -> item[1] }]
                        }
                        .filter { item -> !item[1].isEmpty() }
                        .ifEmpty { error(all_dropped_message) }
        merged = MERGE_OUTPUTS(
            consensus, ref_genome, fai_index, dropped_samples
        )
        EXPORT_OUTPUT(
            merged.map { bcf -> tuple('consensus', bcf) },
            ref_genome
        )
    } else {
        consensus_bcf = GENERATE_CONSENSUS.out.consensus_bcf
            .combine(dropped_report)
            .map { sample_id, bcf, _report -> tuple(sample_id, bcf) }
            .ifEmpty { error(all_dropped_message) }
        EXPORT_OUTPUT(
            consensus_bcf,
            ref_genome
        )
    }
}

