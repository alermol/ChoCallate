#!/usr/bin/env nextflow

include { handleHelpVersion } from './functions/help_version.nf'

handleHelpVersion()

def available_callers         = 'bcftools,gatk,freebayes,snver,vardict'
def diploid_callers           = 'bcftools,gatk,freebayes,snver,vardict'
def polyploid_callers         = 'gatk,freebayes,snver'
def default_diploid_callers   = 'bcftools,gatk,freebayes,snver,vardict'
def default_polyploid_callers = 'gatk,freebayes,snver'
def min_callers_count         = 3

include { 
    getAvailableCallersCount; 
    getConsensusThreshold; 
    normalizeCallerNames;
    allEffectiveCallersInAvailable; 
    effectiveCallersAtLeastThree; 
    allEffectiveCallersDiploidSuitable; 
    allEffectiveCallersPolyploidSuitable;
    validateAllParameters;
    validateTSVFile;
    validateInputFormat;
    validateReferenceGenome;
    validateBowtie2Index;
    validatePloidy;
    validateReadsSource;
    validateConsensusType;
    validateLogLevel;
    validateLogFormat;
    validateWindowSize;
    validateCPUParameter;
    validateForkParameter;
    validateQualityParameter;
    validateNumericParameter
} from './functions/utils.nf'

include { 
    initLogging; 
    logInfo; 
    logWarn; 
    logError; 
    logDebug; 
    logProcessStart; 
    logProcessComplete; 
    logProcessError; 
    logWorkflowStart; 
    logWorkflowComplete; 
    logWorkflowError;
    logPerformance;
    logResourceUsage;
    logFileOperation;
    logValidation
} from './functions/logging.nf'

include { 
    CREATE_FAI_INDEX 
} from './functions/create_fai_index.nf'
include { 
    CREATE_SEQ_DICT 
} from './functions/create_seq_dict.nf'
include { 
    PREPARE_BAM 
} from './functions/prepare_bam.nf'
include { 
    COVERAGE_GENERATION 
} from './functions/coverage_generation.nf'
include { 
    GENERATE_ZERO_BCF 
} from './functions/generate_zero_bcf.nf'
include { 
    CALLING 
} from './functions/calling.nf'
include { 
    GENERATE_CONSENSUS 
} from './functions/generate_consensus.nf'
include { 
    MERGE_BCFS 
} from './functions/merge_bcfs.nf'
include { 
    CLEANUP_SAMPLE_TEMP 
} from './functions/cleanup_sample_temp.nf'

initLogging()

logInfo("Cleanup configuration loaded", [
    enable_sample_cleanup: params.enable_sample_cleanup,
    cleanup_intermediate_bam: params.cleanup_intermediate_bam,
    cleanup_intermediate_bcf: params.cleanup_intermediate_bcf,
    cleanup_intermediate_subfolders: params.cleanup_intermediate_subfolders,
    cleanup_input_symlinks: params.cleanup_input_symlinks,
    debug_mode: params.debug,
    note: "Work directory always persists. In debug mode, all intermediate files are preserved. In production mode, sample-level cleanup is performed."
])

logInfo("Starting input validation and parameter sanity checks")

def validationResult = validateAllParameters(params)

if (!validationResult.valid) {
    logError("Input validation failed", [
        validation: "comprehensive_parameter_validation",
        error_count: validationResult.errors.size(),
        errors: validationResult.errors
    ])
    
    validationResult.errors.each { error ->
        logError("Validation error", [error: error])
    }
    
    error "Input validation failed with ${validationResult.errors.size()} error(s). Please fix the issues above and try again."
}

if (validationResult.warnings.size() > 0) {
    logWarn("Parameter warnings detected", [
        validation: "parameter_warnings",
        warning_count: validationResult.warnings.size(),
        warnings: validationResult.warnings
    ])
    
    validationResult.warnings.each { warning ->
        logWarn("Parameter warning", [warning: warning])
    }
}

logInfo("Input validation completed successfully", [
    validation: "comprehensive_parameter_validation",
    parameters_validated: true,
    warnings_count: validationResult.warnings.size()
])

def outdir = new File(params.outdir)
if (!outdir.exists()) {
    try {
        outdir.mkdirs()
        logInfo("Output directory created", [
            validation: "output_directory_creation",
            directory: params.outdir
        ])
    } catch (Exception e) {
        logError("Failed to create output directory", [
            validation: "output_directory_creation",
            directory: params.outdir,
            error: e.message
        ])
        error "Failed to create output directory: ${params.outdir}. Error: ${e.message}"
    }
} else if (!outdir.canWrite()) {
    logError("Output directory is not writable", [
        validation: "output_directory_permissions",
        directory: params.outdir
    ])
    error "Output directory is not writable: ${params.outdir}"
} else {
    logInfo("Output directory validation passed", [
        validation: "output_directory_permissions",
        directory: params.outdir,
        writable: true
    ])
}

def samplesValidation = validateTSVFile(params.samples_tsv, params.input_format)
if (samplesValidation.valid) {
    logInfo("Samples TSV validation passed", [
        validation: "samples_tsv_format",
        file: params.samples_tsv,
        sample_count: samplesValidation.sampleCount
    ])
}

def refGenomeValidation = validateReferenceGenome(params.reference_genome)
if (refGenomeValidation.valid) {
    logInfo("Reference genome validation passed", [
        validation: "reference_genome_format",
        file: params.reference_genome
    ])
}

if (params.input_format == 'fastq') {
    def refIndexValidation = validateBowtie2Index(params.reference_index)
    if (refIndexValidation.valid) {
        logInfo("Bowtie2 index validation passed", [
            validation: "bowtie2_index_format",
            index_directory: params.reference_index,
            index_files_count: refIndexValidation.indexFiles.size()
        ])
    }
}

def effective_callers
if (params.effective_callers == "-") {
    if (params.ploidy == 2) {
        effective_callers = default_diploid_callers
        logInfo("Using default diploid callers", [callers: effective_callers])
    } else if (params.ploidy > 2) {
        effective_callers = default_polyploid_callers
        logInfo("Using default polyploid callers", [callers: effective_callers])
    } else {
        logError("Invalid ploidy value", [ploidy: params.ploidy, error: "Ploidy must be 2 or greater"])
        error "Invalid ploidy value: ${params.ploidy}. Ploidy must be 2 or greater."
    }
} else {
    effective_callers = normalizeCallerNames(params.effective_callers)
    logInfo("Using user-defined effective callers", [callers: effective_callers, original: params.effective_callers])
}

if (!allEffectiveCallersInAvailable(effective_callers, available_callers)) {
    logError("Validation failed", [validation: "effective_callers_availability", callers: effective_callers, available: available_callers])
    exit 1
} else {
    logInfo("Validation passed", [validation: "effective_callers_availability", callers: effective_callers])
}

if (effectiveCallersAtLeastThree(effective_callers) < min_callers_count) {
    logError("Validation failed", [
        validation: "minimum_callers_count", 
        required: min_callers_count, 
        provided: getAvailableCallersCount(effective_callers), 
        callers: effective_callers
    ])
    exit 1
}

if (params.ploidy == 2) {
    if (allEffectiveCallersDiploidSuitable(effective_callers, diploid_callers)) {
        logInfo("Validation passed", [validation: "diploid_callers_suitability", callers: effective_callers])
    } else {
        logError("Validation failed", [validation: "diploid_callers_suitability", callers: effective_callers, suitable: diploid_callers])
        exit 1
    }
} else if (params.ploidy > 2) {
    if (allEffectiveCallersPolyploidSuitable(effective_callers, polyploid_callers)) {
        logInfo("Validation passed", [validation: "polyploid_callers_suitability", callers: effective_callers])
    } else {
        logError("Validation failed", [validation: "polyploid_callers_suitability", callers: effective_callers, suitable: polyploid_callers])
        exit 1
    }
}

if (params.debug) {
    logInfo("Debug mode enabled")
} else {
    logInfo("Debug mode disabled")
}

workflow {
    logWorkflowStart()
    
    if (params.test_run) {
        logInfo("Test run mode enabled")
        Channel
            .fromPath(params.samples_tsv)
            .splitCsv(header: false, sep: '\t', limit: params.test_run_limit)
            .map{row -> tuple(row[0], row[1], row[2], row[3])}
            .set{sample_run_ch}
    } else {
        Channel
            .fromPath(params.samples_tsv)
            .splitCsv(header: false, sep: '\t')
            .map{row -> tuple(row[0], row[1], row[2], row[3])}
            .set{sample_run_ch}
    }
    
    logInfo("Sample channel created", [samples_count: "processing"])

    ref_index = params.input_format == 'fastq' ? file(params.reference_index) : 'NA'
    ref_genome = file(params.reference_genome)

    CREATE_FAI_INDEX(ref_genome)
    CREATE_SEQ_DICT(CREATE_FAI_INDEX.out.fai_index.collect{it[0]})

    PREPARE_BAM(sample_run_ch,
                CREATE_SEQ_DICT.out.gen_dict,
                CREATE_FAI_INDEX.out.fai_index,
                ref_index)

    COVERAGE_GENERATION(PREPARE_BAM.out.bam)

    GENERATE_ZERO_BCF(PREPARE_BAM.out.bam,
                      CREATE_FAI_INDEX.out.fai_index,
                      COVERAGE_GENERATION.out.coverage)

    def parallel_cpus = getAvailableCallersCount(effective_callers)
    CALLING(PREPARE_BAM.out.bam,
            CREATE_FAI_INDEX.out.fai_index,
            COVERAGE_GENERATION.out.coverage,
            params.ploidy,
            CREATE_SEQ_DICT.out.gen_dict,
            effective_callers,
            parallel_cpus,
            params.bcftools_mpileup_extra_args,
            params.bcftools_call_extra_args,
            params.freebayes_extra_args,
            params.gatk4_extra_args,
            params.vardict_extra_args,
            params.snver_extra_args)

    def cons_threshold = getConsensusThreshold(params.cons_type, available_callers)
    GENERATE_CONSENSUS(CALLING.out.snps_vcf,
                       CALLING.out.indels_vcf,
                       CREATE_FAI_INDEX.out.fai_index.collect{it[1]},
                       GENERATE_ZERO_BCF.out.zero_bcf,
                       cons_threshold)

    MERGE_BCFS(GENERATE_CONSENSUS.out.final_snps.map{it[1]}.collect(), 
               GENERATE_CONSENSUS.out.final_indels.map{it[1]}.collect())

    if (params.enable_sample_cleanup && !params.debug) {
        logInfo("Sample-level cleanup enabled - cleaning intermediate files", [
            action: "sample_cleanup_enabled",
            debug_mode: false,
            intermediate_bam_cleanup: params.cleanup_intermediate_bam,
            intermediate_vcf_cleanup: params.cleanup_intermediate_bcf
        ])
        
        CLEANUP_SAMPLE_TEMP(GENERATE_CONSENSUS.out.final_snps,
                            GENERATE_CONSENSUS.out.final_indels,
                            PREPARE_BAM.out.bam,
                            COVERAGE_GENERATION.out.coverage,
                            GENERATE_ZERO_BCF.out.zero_bcf,
                            CALLING.out.snps_vcf,
                            CALLING.out.indels_vcf,
                            MERGE_BCFS.out.merged,
                            params.cleanup_intermediate_bam,
                            params.cleanup_intermediate_bcf)
    } else if (params.enable_sample_cleanup && params.debug) {
        logInfo("Sample-level cleanup skipped in debug mode - preserving all intermediate files", [
            action: "sample_cleanup_skipped",
            debug_mode: true,
            reason: "Debug mode preserves intermediate files for analysis",
            intermediate_bam_cleanup: params.cleanup_intermediate_bam,
            intermediate_vcf_cleanup: params.cleanup_intermediate_bcf
        ])
    } else if (!params.enable_sample_cleanup) {
        logInfo("Sample-level cleanup disabled by configuration", [
            action: "sample_cleanup_disabled",
            reason: "enable_sample_cleanup = false"
        ])
    }
}

workflow.onComplete {
    logWorkflowComplete()
    
    if (params.enable_sample_cleanup && !params.debug) {
        logInfo("Sample-level cleanup was executed during pipeline execution", [
            action: "cleanup_summary",
            sample_cleanup_enabled: true,
            debug_mode: false,
            status: "intermediate_files_cleaned_via_sample_cleanup"
        ])
    } else if (params.enable_sample_cleanup && params.debug) {
        logInfo("Sample-level cleanup was skipped due to debug mode", [
            action: "cleanup_summary",
            sample_cleanup_enabled: true,
            debug_mode: true,
            status: "intermediate_files_preserved_for_debugging",
            reason: "Debug mode preserves all intermediate files for analysis and troubleshooting"
        ])
    } else if (!params.enable_sample_cleanup) {
        logInfo("Sample-level cleanup was disabled by configuration", [
            action: "cleanup_summary",
            sample_cleanup_enabled: false,
            status: "no_sample_cleanup_performed",
            reason: "enable_sample_cleanup = false"
        ])
    }
    
    def workDir = workflow.workDir ? file(workflow.workDir) : null
    if (workDir?.exists()) {
        if (params.debug) {
            logInfo("Debug mode cleanup - preserving entire work directory for debugging", [
                action: "debug_cleanup_preserve",
                workDir: workDir.toString(),
                reason: "Debug mode preserves all files for analysis and troubleshooting"
            ])
        } else {
            logInfo("Production mode cleanup - work directory preserved, sample-level cleanup completed", [
                action: "production_cleanup_complete",
                workDir: workDir.toString(),
                note: "Work directory preserved for analysis. Sample-level cleanup was performed during pipeline execution to clean intermediate files."
            ])
        }
        
        logInfo("Workflow cleanup completed - work directory preserved", [
            action: "cleanup_complete", 
            workDir: workDir.toString(),
            result: "work_directory_preserved",
            note: "All intermediate files are either cleaned (production mode) or preserved (debug mode)"
        ])
    }
}

workflow.onError {
    logWorkflowError(workflow.errorMessage)
}




