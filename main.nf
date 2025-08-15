#!/usr/bin/env nextflow

def available_callers         = 'bcftools,gatk,freebayes,snver,vardict'
def diploid_callers           = 'bcftools,gatk,freebayes,snver,vardict'
def polyploid_callers         = 'gatk,freebayes,snver'
def default_diploid_callers   = 'bcftools,gatk,freebayes,snver,vardict'
def default_polyploid_callers = 'gatk,freebayes,snver'
def min_callers_count         = 3


include { 
    getAvailableCallersCount; 
    getConsensusThreshold; 
    allEffectiveCallersInAvailable; 
    effectiveCallersAtLeastThree; 
    allEffectiveCallersDiploidSuitable; 
    allEffectiveCallersPolyploidSuitable;
    validateAllParameters;
    validateTSVFile;
    validateReferenceGenome;
    validateBowtie2Index;
    validatePloidy;
    validateReadsType;
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

// Include logging module
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

// Initialize logging
initLogging()

// Log cleanup configuration
    logInfo("Cleanup configuration loaded", [
        enable_sample_cleanup: params.enable_sample_cleanup,
        cleanup_intermediate_bam: params.cleanup_intermediate_bam,
        cleanup_intermediate_vcf: params.cleanup_intermediate_vcf,
        note: "Work directory always persists. Sample-level cleanup performed in production mode to clean intermediate files."
    ])

// =============================================================================
// INPUT VALIDATION AND PARAMETER SANITY CHECKS
// =============================================================================

logInfo("Starting input validation and parameter sanity checks")

// Validate all parameters comprehensively
def validationResult = validateAllParameters(params)

if (!validationResult.valid) {
    logError("Input validation failed", [
        validation: "comprehensive_parameter_validation",
        error_count: validationResult.errors.size(),
        errors: validationResult.errors
    ])
    
    // Log each error individually for clarity
    validationResult.errors.each { error ->
        logError("Validation error", [error: error])
    }
    
    error "Input validation failed with ${validationResult.errors.size()} error(s). Please fix the issues above and try again."
}

// Log warnings if any
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

// Validate output directory
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

// Validate specific file types individually for detailed logging
def samplesValidation = validateTSVFile(params.samples_tsv)
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

def refIndexValidation = validateBowtie2Index(params.reference_index)
if (refIndexValidation.valid) {
    logInfo("Bowtie2 index validation passed", [
        validation: "bowtie2_index_format",
        index_directory: params.reference_index,
        index_files_count: refIndexValidation.indexFiles.size()
    ])
}

// =============================================================================
// EXISTING CALLER VALIDATION (ENHANCED WITH LOGGING)
// =============================================================================

// Handle effective_callers parameter
// If params.effective_callers has default value "-", then set of caller in default_diploid_callers or default_polyploid_callers will be taken as default for diploid and polyploid calling
// If params.effective_callers have user-defined value, then it should be used
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
    effective_callers = params.effective_callers
    logInfo("Using user-defined effective callers", [callers: effective_callers])
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


// Main workflow definition
workflow {
    logWorkflowStart()
    
    // Create channel from samples TSV file and parse into tuples
    Channel
        .fromPath(params.samples_tsv)
        .splitCsv(header: false, sep: '\t')
        .map{row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3]))}
        .set{sample_run_ch}
    
    logInfo("Sample channel created", [samples_count: "processing"])

    // Define reference files
    ref_index = file(params.reference_index)
    ref_genome = file(params.reference_genome)

    // Create reference genome index files
    CREATE_FAI_INDEX(ref_genome)
    CREATE_SEQ_DICT(ref_genome)

    // Map reads to reference genome
    PREPARE_BAM(sample_run_ch,
                CREATE_SEQ_DICT.out.gen_dict,
                CREATE_FAI_INDEX.out.fai_index,
                ref_index)

    // Generate coverage information from BAM files
    COVERAGE_GENERATION(PREPARE_BAM.out.bam)

    GENERATE_ZERO_VCF(PREPARE_BAM.out.bam,
                      CREATE_FAI_INDEX.out.fai_index,
                      COVERAGE_GENERATION.out.coverage)

    // Call variants
    CALLING(PREPARE_BAM.out.bam,
            CREATE_FAI_INDEX.out.fai_index,
            COVERAGE_GENERATION.out.coverage,
            params.ploidy,
            CREATE_SEQ_DICT.out.gen_dict)

    // Generate final consensus VCF
    GENERATE_CONSENSUS(CALLING.out.snps_vcf,
                       CALLING.out.indels_vcf,
                       CREATE_FAI_INDEX.out.fai_index.map{it[1]},
                       GENERATE_ZERO_VCF.out.zero_vcf)
    
    // Cleanup temporary files for each sample after processing is complete (if enabled and NOT in debug mode)
    if (params.enable_sample_cleanup && !params.debug) {
        logInfo("Sample-level cleanup enabled - cleaning intermediate files", [
            action: "sample_cleanup_enabled",
            debug_mode: false,
            intermediate_bam_cleanup: params.cleanup_intermediate_bam,
            intermediate_vcf_cleanup: params.cleanup_intermediate_vcf
        ])
        
        CLEANUP_SAMPLE_TEMP(GENERATE_CONSENSUS.out.final_snps,
                            GENERATE_CONSENSUS.out.final_indels,
                            PREPARE_BAM.out.bam,
                            COVERAGE_GENERATION.out.coverage,
                            GENERATE_ZERO_VCF.out.zero_vcf,
                            CALLING.out.snps_vcf,
                            CALLING.out.indels_vcf,
                            params.cleanup_intermediate_bam,
                            params.cleanup_intermediate_vcf)
    } else if (params.enable_sample_cleanup && params.debug) {
        logInfo("Sample-level cleanup skipped in debug mode - preserving all intermediate files", [
            action: "sample_cleanup_skipped",
            debug_mode: true,
            reason: "Debug mode preserves intermediate files for analysis",
            intermediate_bam_cleanup: params.cleanup_intermediate_bam,
            intermediate_vcf_cleanup: params.cleanup_intermediate_vcf
        ])
    } else if (!params.enable_sample_cleanup) {
        logInfo("Sample-level cleanup disabled by configuration", [
            action: "sample_cleanup_disabled",
            reason: "enable_sample_cleanup = false"
        ])
    }
}



// Cleanup temporary files after workflow completion
workflow.onComplete {
    logWorkflowComplete()
    
    // Log cleanup summary based on configuration and debug mode
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
    
    // Workflow-level cleanup
    def workDir = workflow.workDir ? file(workflow.workDir) : null
    if (workDir?.exists()) {
        if (params.debug) {
            // In debug mode, preserve entire work directory - no cleanup performed
            logInfo("Debug mode cleanup - preserving entire work directory for debugging", [
                action: "debug_cleanup_preserve",
                workDir: workDir.toString(),
                reason: "Debug mode preserves all files for analysis and troubleshooting"
            ])
        } else {
            // In production mode, preserve work directory but sample-level cleanup was already performed
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

// Process to create FASTA index file (for FreeBayes)
process CREATE_FAI_INDEX {
    maxForks 1
    cpus 1

    input:
    path(ref_genome)

    output:
    tuple path("${ref_genome}"), path("${ref_genome}.fai"), emit: fai_index

    script:
    def startTime = System.currentTimeMillis()
    def refName = ref_genome.getName()
    
    """
    # Log process start
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Process started - Creating FASTA index"
    
    samtools faidx --threads ${task.cpus} ${ref_genome}
    
    # Log process completion and performance
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Process completed - FASTA index created"
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Performance: completed successfully"
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
    def startTime = System.currentTimeMillis()
    def refName = ref_genome.getName()
    
    """
    # Log process start
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Process started - Creating sequence dictionary"
    
    gatk CreateSequenceDictionary -R ${ref_genome}

    # Remove input file after completion
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Removing input file: ${ref_genome}"
    rm -f ${ref_genome} 2>/dev/null || true
    
    # Log process completion and performance
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Process completed - Sequence dictionary created"
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Performance: completed successfully"
    """
}

// Process to map reads to reference genome
process PREPARE_BAM {
    maxForks params.bowtie2_forks
    cpus params.bowtie2_cpu

    tag "${sample_id}"

    input:
    tuple val(sample_id), val(read1), val(read2), val(read3)
    path(genome_dictionary)
    tuple path(ref_genome), path(genome_fai)
    val(ref_index)

    output:
    tuple path("${sample_id}.bam"), path("${sample_id}.bam.csi"), emit: bam

    script:
    def startTime = System.currentTimeMillis()
    def readsType = params.reads_type
    
    if ( params.reads_type == 'pe' )
        """
        # Log process start
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process started - Mapping reads (${readsType})"
        
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -1 ${read1} -2 ${read2} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.primary.bam

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Running LeftAlignIndels on primary BAM"
        gatk LeftAlignIndels -I ${sample_id}.primary.bam -O ${sample_id}.bam -R ${ref_genome} -OBI false

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Removing intermediate primary BAM"
        rm -f ${sample_id}.primary.bam

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Indexing final BAM"
        samtools index --csi --threads ${task.cpus} ${sample_id}.bam

        # Remove input files after completion
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Removing input files: reads and reference"
        rm -f ${ref_genome} ${genome_dictionary} ${genome_fai} 2>/dev/null || true

            
        # Log process completion and performance
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process completed - BAM file created"
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Performance: completed successfully"
        """
    else if ( params.reads_type == 'se' )
        """
        # Log process start
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process started - Mapping reads (${readsType})"
        
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -U ${read1} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.primary.bam

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Running LeftAlignIndels on primary BAM"
        gatk LeftAlignIndels -I ${sample_id}.primary.bam -O ${sample_id}.bam -R ${ref_genome} -OBI false

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Removing intermediate primary BAM"
        rm -f ${sample_id}.primary.bam

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Indexing final BAM"
        samtools index --csi --threads ${task.cpus} ${sample_id}.bam

        # Remove input files after completion
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Removing input files: reads and reference"
        rm -f ${ref_genome} ${genome_dictionary} ${genome_fai} 2>/dev/null || true
            
        # Log process completion and performance
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process completed - BAM file created"
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Performance: completed successfully"
        """
    else if ( params.reads_type == 'mx' )
        """
        # Log process start
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process started - Mapping reads (${readsType})"
        
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -1 ${read1} -2 ${read2} -U ${read3} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.primary.bam

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Running LeftAlignIndels on primary BAM"
        gatk LeftAlignIndels -I ${sample_id}.primary.bam -O ${sample_id}.bam -R ${ref_genome} -OBI false

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Removing intermediate primary BAM"
        rm -f ${sample_id}.primary.bam

        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Indexing final BAM"
        samtools index --csi --threads ${task.cpus} ${sample_id}.bam

        # Remove input files after completion
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Removing input files: reads and reference"
        rm -f ${ref_genome} ${genome_dictionary} ${genome_fai} 2>/dev/null || true
            
        # Log process completion and performance
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Process completed - BAM file created"
        echo "[\$(date -Iseconds)] [INFO] [BOWTIE2_MAPPING] [${sample_id}] Performance: completed successfully"
        """
    else
        error 'Invalid reads type: ${params.reads_type}. Available types: se, pe, mx'
}

// Process to generate coverage information
process COVERAGE_GENERATION {
    cpus 1
    maxForks 1

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)

    output:
    path("${bam.baseName}.bed"), emit: coverage

    script:
    def startTime = System.currentTimeMillis()
    def bamName = bam.getName()
    
    """
    # Log process start
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bamName}] Process started - Generating coverage information"
    
    samtools depth -J --threads ${task.cpus} ${bam} | \
        awk '\$3 >= ${params.min_coverage} {print \$1,\$2-1,\$2}' | \
        bedops --merge - > ${bam.baseName}.bed
    
    # Remove input files after completion
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bamName}] Removing input files: BAM and BAM index"
    rm -f ${bam} ${bam_index}
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bamName}] Input files removed"

    # Log process completion and performance
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bamName}] Process completed - Coverage BED file created"
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bamName}] Performance: completed successfully"
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
    def startTime = System.currentTimeMillis()
    def bamName = bam.getName()
    String genotype = (["0"] * params.ploidy).join("/")
    
    """
    # Log process start
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_VCF] [${bamName}] Process started - Generating zero VCF"
    
    bcftools mpileup -Ov --count-orphans --fasta-ref ${ref_genome} --threads ${task.cpus} --max-depth 1 \
        --min-BQ ${params.min_base_quality} --regions-file ${coverage_bed} ${bam} | \
        awk -v OFS='\\t' -v gen=${genotype} '{if(\$0 !~ /#/) print \$1,\$2,\$3,\$4,".","100",".",".","GT",gen; else print \$0}' | \
        awk -v OFS='\\t' '{if(length(\$4) == 1 || \$0 ~ /#/) print \$0}' | bgzip  > zero.vcf.gz

    # Remove input files after completion
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_VCF] [${bamName}] Removing input files: BAM, BAM index, reference genome, reference genome FAI, and coverage BED"
    rm -f ${bam} ${bam_index} ${ref_genome} ${ref_genome_fai} ${coverage_bed}
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_VCF] [${bamName}] Input files removed"
    
    # Log process completion and performance
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_VCF] [${bamName}] Process completed - Zero VCF created"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_VCF] [${bamName}] Performance: completed successfully"
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
    def startTime = System.currentTimeMillis()
    def bamName = bam.getName()
    def parallel_cpus = getAvailableCallersCount(effective_callers)
    
    """
    # Log process start
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Process started - Running variant callers"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Using callers: ${effective_callers}"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Parallel CPUs: ${parallel_cpus}"
    
    touch callers_commands.sh

    if [[ ",${effective_callers}," == *"bcftools"* ]]; then
        echo "bash ${projectDir}/bin/bcftools_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.bcftools_cpu}" >> callers_commands.sh
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Added bcftools caller command"
    fi

    if [[ ",${effective_callers}," == *"freebayes"* ]]; then
        echo "bash ${projectDir}/bin/freebayes_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.reads_source} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Added freebayes caller command"
    fi

    if [[ ",${effective_callers}," == *"gatk"* ]]; then
        echo "bash ${projectDir}/bin/gatk4_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Added GATK4 caller command"
    fi

    if [[ ",${effective_callers}," == *"vardict"* ]]; then
        echo "bash ${projectDir}/bin/vardict_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.vardict_cpu}" >> callers_commands.sh
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Added VarDict caller command"
    fi

    if [[ ",${effective_callers}," == *"snver"* ]]; then
        echo "bash ${projectDir}/bin/snver_caller.sh ${bam} ${coverage_bed} ${ref_genome} ${ploidy} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Added SNVer caller command"
    fi

    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Executing ${parallel_cpus} callers in parallel"
    parallel -j ${parallel_cpus} '{}' :::: callers_commands.sh
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Removing callers_commands.sh and all input/intermediate files"
    # Remove callers_commands.sh
    rm -f callers_commands.sh 2>/dev/null || true

    # Remove all input and intermediate files
    rm -f "${bam}" 2>/dev/null || true
    rm -f "${bam_index}" 2>/dev/null || true
    rm -f "${coverage_bed}" 2>/dev/null || true
    rm -f "${ref_genome}" 2>/dev/null || true
    rm -f "${ref_genome_fai}" 2>/dev/null || true
    rm -f "${ref_genome_dict}" 2>/dev/null || true
    rm -f "${bam.baseName}.raw.vcf" 2>/dev/null || true
    rm -f "${bam.baseName}.indel.raw.vcf" 2>/dev/null || true
    rm -f "${bam.baseName}.failed.log" 2>/dev/null || true
    rm -f "${bam.baseName}.gatk1.vcf.gz.csi" 2>/dev/null || true

    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] callers_commands.sh, input files, and intermediate files removed"

    # Log process completion and performance
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Process completed - Variant calling finished"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bamName}] Performance: completed successfully"
    """
}


// Process to generate final consensus VCF
process GENERATE_CONSENSUS {
    maxForks params.cons_forks
    cpus params.cons_cpus

    tag "${sample}"

    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.snps.vcf.gz'
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.indels.vcf.gz'

    input:
    tuple val(sample), path("${sample}.snps_*.vcf.gz", arity: '3..*')
    tuple val(sample), path("${sample}.indels_*.vcf.gz", arity: '3..*')
    path(ref_genome_fai)
    path(zero_vcf)

    output:
    tuple val("${sample}"), path("${sample}.snps.vcf.gz"), emit: final_snps
    tuple val("${sample}"), path("${sample}.indels.vcf.gz"), emit: final_indels

    script:
    def startTime = System.currentTimeMillis()
    def cons_threshold = getConsensusThreshold(params.cons_type, available_callers)
    
    """
    # Log process start
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Process started - Generating consensus VCFs"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Consensus threshold: ${cons_threshold}"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Window size: ${params.win_size}"
    
    awk -v OFS='\t' '{print \$1,"0",\$2}' ${ref_genome_fai} > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Processing SNPs..."
    for i in ${sample}.snps_*; do tabix -C \${i}; done

    tabix -C zero.vcf.gz

    parallel -j ${task.cpus} 'consensus_generation.sh {1} {#} ${sample} "snps" ${cons_threshold}' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.snps.vcf.gz

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] SNPs consensus VCF created"
    rm -r all_chrs/*

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Processing indels..."
    for i in ${sample}.indels_*; do tabix -C \${i}; done

    parallel -j ${task.cpus} 'consensus_generation.sh {1} {#} ${sample} "indels" ${cons_threshold}' :::: genome_intervals.bed

    find all_chrs/ -name '*.vcf.gz' -type f > vcf_files.txt

    bcftools concat --naive-force -Oz --file-list vcf_files.txt | \
        bcftools reheader --threads ${task.cpus} -f ${ref_genome_fai} | \
        bcftools sort -Oz -o ${sample}.indels.vcf.gz
        
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indels consensus VCF created"
    
    # Cleanup intermediate and input files (including those not deleted by consensus_generation.sh)
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Cleaning up intermediate and input files"

    # Remove BED and VCF list files
    rm -f genome.bed genome_intervals.bed vcf_files.txt 2>/dev/null || true

    # Remove all_chrs directory if not already removed by consensus_generation.sh
    rm -rf all_chrs 2>/dev/null || true

    # Remove input VCFs and their indexes for SNPs and indels
    rm -f ${sample}.snps_* 2>/dev/null || true
    rm -f ${sample}.indels_* 2>/dev/null || true

    # Remove zero.vcf.gz and its index if present
    rm -f zero.vcf.gz zero.vcf.gz.csi 2>/dev/null || true

    # Remove fai index file
    rm -f ${ref_genome_fai} 2>/dev/null || true

    # Log process completion and performance
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Process completed - Consensus VCFs generated"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Performance: completed successfully"
    """
}

// Process to cleanup all temporary files for a specific sample
// Each process runs in its own unique directory, avoiding file name collisions
// IMPORTANT: Files specified as input parameters are preserved and NOT deleted
// Symbolic links are followed and the actual target files are deleted
process CLEANUP_SAMPLE_TEMP {
    tag "${sample}"
    
    input:
    tuple val(sample), path(output_vcf_snps)
    tuple val(sample), path(output_vcf_indels)
    tuple path(bam), path(bam_index)
    path(coverage_bed)
    path(zero_vcf)
    tuple val(sample), path("${sample}.snps_*.vcf.gz", arity: '3..*')
    tuple val(sample), path("${sample}.indels_*.vcf.gz", arity: '3..*')
    val(cleanup_intermediate_bam)
    val(cleanup_intermediate_vcf)

    // output:
    // tuple val("${sample_id}"), val("cleaned"), emit: cleanup_status
    
    script:
    """
    log_msg() {
        local level="\$1"
        local msg="\$2"
        echo "[\$(date -Iseconds)] [${sample}] [CLEANUP_SAMPLE_TEMP] [\${level}] \${msg}"
    }

    log_msg "INFO" "Starting sample cleanup process (symlink-aware)"

    # Function to remove a file or symlink, following symlinks to the real file
    remove_file_follow_symlink() {
        local f="\$1"
        if [ -L "\${f}" ]; then
            # If it's a symlink, resolve and remove the target
            local target
            target=\$(realpath "\${f}")
            if [ -e "\$target" ]; then
                rm -f "\$target" && log_msg "INFO" "Removed symlink target: \$target"
            fi
            rm -f "\$f" && log_msg "INFO" "Removed symlink: \$f"
        elif [ -e "\$f" ]; then
            rm -f "\$f" && log_msg "INFO" "Removed file: \$f"
        fi
    }

    # Remove SNP and indel input VCFs and their indexes (all matching files)
    for vcf in ${output_vcf_snps}; do
        [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        [ -e "\${vcf}.csi" ] && remove_file_follow_symlink "\${vcf}.csi"
    done

    for vcf in ${output_vcf_indels}; do
        [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        [ -e "\${vcf}.csi" ] && remove_file_follow_symlink "\${vcf}.csi"
    done

    # Remove BAM and BAM index if cleanup_intermediate_bam is true
    if [ "${cleanup_intermediate_bam}" = "true" ]; then
        [ -e "${bam}" ] && remove_file_follow_symlink "${bam}"
        [ -e "${bam_index}" ] && remove_file_follow_symlink "${bam_index}"
        log_msg "INFO" "Removed BAM and BAM index files"
    else
        log_msg "INFO" "Preserved BAM and BAM index files (cleanup_intermediate_bam = false)"
    fi

    # Remove coverage BED file
    [ -e "${coverage_bed}" ] && remove_file_follow_symlink "${coverage_bed}"

    # Remove zero VCF
    [ -e "${zero_vcf}" ] && remove_file_follow_symlink "${zero_vcf}"

    # Remove consensus VCFs if cleanup_intermediate_vcf is true
    if [ "${cleanup_intermediate_vcf}" = "true" ]; then
        # Remove all consensus VCFs for each caller (using globbing, since effective_callers may not be a number)
        for vcf in "${sample}.snps_"*.vcf.gz "${sample}.indels_"*.vcf.gz; do
            [ -e "\$vcf" ] && remove_file_follow_symlink "\$vcf"
        done
        log_msg "INFO" "Removed consensus VCFs"
    else
        log_msg "INFO" "Preserved consensus VCFs (cleanup_intermediate_vcf = false)"
    fi

    # Remove all symbolic links in the current directory
    for symlink in ./*; do
        if [ -L "\$symlink" ]; then
            rm -f "\$symlink" && log_msg "INFO" "Removed symlink in current directory: \$symlink"
        fi
    done

    log_msg "INFO" "Sample cleanup completed"
    """
}
