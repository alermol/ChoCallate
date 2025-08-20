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
    validateInputFormat;
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

def samplesValidation = validateTSVFile(params.samples_tsv, params.input_format, params.reads_type)
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

workflow {
    logWorkflowStart()
    
    Channel
        .fromPath(params.samples_tsv)
        .splitCsv(header: false, sep: '\t')
        .map{row -> tuple(row[0], row[1], row[2], row[3])}
        .set{sample_run_ch}
    
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

    CALLING(PREPARE_BAM.out.bam,
            CREATE_FAI_INDEX.out.fai_index,
            COVERAGE_GENERATION.out.coverage,
            params.ploidy,
            CREATE_SEQ_DICT.out.gen_dict)

    GENERATE_CONSENSUS(CALLING.out.snps_vcf,
                       CALLING.out.indels_vcf,
                       CREATE_FAI_INDEX.out.fai_index.collect{it[1]},
                       GENERATE_ZERO_BCF.out.zero_bcf)
    
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

process CREATE_FAI_INDEX {
    maxForks 1
    cpus 1

    input:
    path(ref_genome)

    output:
    tuple path("reference.fasta"), path("reference.fasta.fai"), emit: fai_index

    script:
    def refName = ref_genome.getName()
    def isGzipped = refName.endsWith('.gz')
    
    """
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Process started - Creating FASTA index"
    
    if [[ "$isGzipped" == "true" ]]; then
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Detected gzipped reference genome, creating ungzipped version"
        gunzip -c ${ref_genome} > reference.fasta
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Ungzipped version created"
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Creating FASTA index"
        samtools faidx --threads ${task.cpus} reference.fasta
    else
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Reference genome is uncompressed"
        mv ${ref_genome} reference.fasta
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Creating FASTA index"
        samtools faidx --threads ${task.cpus} reference.fasta
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Process completed - FASTA index created"
    """
}

process CREATE_SEQ_DICT {
    maxForks 1

    input:
    path(ref_genome)

    output:
    path("${ref_genome.baseName}.dict"), emit: gen_dict

    script:
    def refName = ref_genome.getName()
    
    """
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Process started - Creating sequence dictionary"
    
    gatk CreateSequenceDictionary -R ${ref_genome}

    rm -f ${ref_genome} 2>/dev/null || true
    
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Process completed - Sequence dictionary created"
    """
}

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
    """
    prepare_bam.sh \
        ${params.input_format} \
        ${params.reads_type} \
        ${sample_id} \
        "${read1}" \
        "${read2}" \
        "${read3}" \
        "${ref_genome}" \
        "${genome_dictionary}" \
        "${genome_fai}" \
        "${ref_index}" \
        ${params.min_map_qual} \
        ${task.cpus} \
        "${params.cleanup_intermediate_subfolders}" \
        "${params.cleanup_input_symlinks}"
    """
}

process COVERAGE_GENERATION {
    cpus 1
    maxForks 1

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)

    output:
    path("${bam.baseName}.bed"), emit: coverage

    script:
    """
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Process started - Generating coverage information"
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Creating intermediate subfolder for coverage generation"
    mkdir -p "${bam.baseName}_coverage_gen"
    cd "${bam.baseName}_coverage_gen"
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Creating symlinks to input files"
    ln -sf "../${bam}" input.bam
    ln -sf "../${bam_index}" input.bam.csi
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Running samtools depth with min coverage ${params.min_coverage}"
    samtools depth -J --threads ${task.cpus} input.bam | \
        awk '\$3 >= ${params.min_coverage} {print \$1,\$2-1,\$2}' | \
        bedops --merge - > ${bam.baseName}.bed
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Moving final output to parent directory"
    mv ${bam.baseName}.bed ../${bam.baseName}.bed
    
    cd ..
    
    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${bam.baseName}_coverage_gen"
        echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${bam}" "${bam_index}" 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [COVERAGE_GENERATION] [${bam.baseName}] Process completed - Coverage BED file created"
    """
}

process GENERATE_ZERO_BCF {
    maxForks params.zero_bcf_forks
    cpus params.zero_bcf_cpu

    tag "${bam.baseName}"

    input:
    tuple path(bam), path(bam_index)
    tuple path(ref_genome), path(ref_genome_fai)
    path(coverage_bed)

    output:
    path("zero.bcf"), emit: zero_bcf

    script:
    String genotype = (["0"] * params.ploidy).join("/")
    
    """
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Process started - Generating zero BCF"
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Creating intermediate subfolder for zero BCF generation"
    mkdir -p "${bam.baseName}_zero_bcf_gen"
    cd "${bam.baseName}_zero_bcf_gen"
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Creating symlinks to input files"
    ln -sf "../${bam}" input.bam
    ln -sf "../${bam_index}" input.bam.csi
    ln -sf "../${ref_genome}" ref_genome.fasta
    ln -sf "../${ref_genome_fai}" ref_genome.fasta.fai
    ln -sf "../${coverage_bed}" coverage.bed

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Writing minimal VCF header for bcftools reheader fallback"
    cat > minimal_header.txt <<EOF
    ##fileformat=VCFv4.3
    ##source=ChoCallateZeroVCF
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
    EOF
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Running bcftools mpileup with min base quality ${params.min_base_quality}"
    bcftools mpileup -Ov --count-orphans --fasta-ref ref_genome.fasta --threads ${task.cpus} --max-depth 1 \
        --min-BQ ${params.min_base_quality} --regions-file coverage.bed input.bam | \
        awk -v OFS='\\t' -v gen=${genotype} '{if(\$0 !~ /#/) print \$1,\$2,\$3,\$4,".","100",".",".","GT",gen; else print \$0}' | \
        awk -v OFS='\\t' '{if(length(\$4) == 1 || \$0 ~ /#/) print \$0}' | \
        bcftools reheader -h minimal_header.txt -f ref_genome.fasta.fai - | \
        bcftools view -Ob -o zero.bcf

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Moving final output to parent directory"
    mv zero.bcf ../zero.bcf
    
    cd ..
    
    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${bam.baseName}_zero_bcf_gen"
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${bam}" "${bam_index}" "${ref_genome}" "${ref_genome_fai}" "${coverage_bed}" 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_ZERO_BCF] [${bam.baseName}] Process completed - Zero BCF created"
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
    tuple val("${bam.baseName}"), path("*.snps_*.bcf"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("*.indels_*.bcf"), emit: indels_vcf

    script:
    def parallel_cpus = getAvailableCallersCount(effective_callers)
    
    """
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Process started - Running variant callers"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Using callers: ${effective_callers}"
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Parallel CPUs: ${parallel_cpus}"
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Creating intermediate subfolder for variant calling"
    mkdir -p "${bam.baseName}_calling"
    cd "${bam.baseName}_calling"
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Creating symlinks to input files"
    ln -sf "../${bam}" input.bam
    ln -sf "../${bam_index}" input.bam.csi
    ln -sf "../${ref_genome}" reference.fasta
    ln -sf "../${ref_genome_fai}" reference.fasta.fai
    ln -sf "../${coverage_bed}" coverage.bed
    ln -sf "../${ref_genome_dict}" reference.dict
    
    touch callers_commands.sh

    if [[ ",${effective_callers}," == *"bcftools"* ]]; then
        echo "bash ${projectDir}/bin/bcftools_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.bcftools_cpu}" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"freebayes"* ]]; then
        echo "bash ${projectDir}/bin/freebayes_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.reads_source} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"gatk"* ]]; then
        echo "bash ${projectDir}/bin/gatk4_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"vardict"* ]]; then
        echo "bash ${projectDir}/bin/vardict_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual} ${params.vardict_cpu}" >> callers_commands.sh
    fi

    if [[ ",${effective_callers}," == *"snver"* ]]; then
        echo "bash ${projectDir}/bin/snver_caller.sh input.bam coverage.bed reference.fasta ${ploidy} ${params.min_base_quality} ${params.min_snp_qual}" >> callers_commands.sh
    fi

    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Executing ${parallel_cpus} callers in parallel"
    parallel -j ${parallel_cpus} '{}' :::: callers_commands.sh
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Moving final outputs to parent directory"
    mv *.snps_*.bcf ../
    mv *.indels_*.bcf ../
    
    cd ..
    
    if [ "${params.cleanup_intermediate_subfolders}" = "true" ]; then
        rm -rf "${bam.baseName}_calling"
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Cleaned up intermediate subfolder"
    fi

    if [ "${params.cleanup_input_symlinks}" = "true" ]; then
        rm -f "${bam}" "${bam_index}" "${ref_genome}" "${ref_genome_fai}" "${coverage_bed}" "${ref_genome_dict}" 2>/dev/null || true
        echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Cleaned up input files"
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [CALLING] [${bam.baseName}] Process completed - Variant calling finished"
    """
}

process GENERATE_CONSENSUS {
    maxForks params.cons_forks
    cpus params.cons_cpus

    tag "${sample}"

    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.snps.bcf'
    publishDir "${params.outdir}/${sample}/", mode: 'copy', pattern: '*.indels.bcf'

    input:
    tuple val(sample), path("${sample}.snps_*.bcf", arity: '3..*')
    tuple val(sample), path("${sample}.indels_*.bcf", arity: '3..*')
    path(ref_genome_fai)
    path(zero_bcf)

    output:
    tuple val("${sample}"), path("${sample}.snps.bcf"), emit: final_snps
    tuple val("${sample}"), path("${sample}.indels.bcf"), emit: final_indels

    script:
    def cons_threshold = getConsensusThreshold(params.cons_type, available_callers)
    
    """
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Process started - Generating consensus BCFs"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Consensus threshold: ${cons_threshold}"
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Window size: ${params.win_size}"
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Creating intermediate subfolder for consensus generation"
    mkdir -p "${sample}_consensus_gen"
    cd "${sample}_consensus_gen"
    
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
    
    awk -v OFS='\\t' '{print \$1,"0",\$2}' ref_genome.fasta.fai > genome.bed
    bedtools makewindows -b genome.bed -w ${params.win_size} > genome_intervals.bed

    mkdir all_chrs

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Processing SNPs..."
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indexing SNP BCF files with bcftools"
    for i in ${sample}.snps_*; do bcftools index \${i}; done

    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Indexing zero BCF with bcftools"
    bcftools index zero.bcf

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
    
    echo "[\$(date -Iseconds)] [INFO] [GENERATE_CONSENSUS] [${sample}] Moving final outputs to parent directory"
    mv ${sample}.snps.bcf ../${sample}.snps.bcf
    mv ${sample}.indels.bcf ../${sample}.indels.bcf
    
    cd ..

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

process CLEANUP_SAMPLE_TEMP {
    tag "${sample}"
    
    input:
    tuple val(sample), path(output_vcf_snps)
    tuple val(sample), path(output_vcf_indels)
    tuple path(bam), path(bam_index)
    path(coverage_bed)
    path(zero_bcf)
    tuple val(sample), path("${sample}.snps_*.bcf", arity: '3..*')
    tuple val(sample), path("${sample}.indels_*.bcf", arity: '3..*')
    val(cleanup_intermediate_bam)
    val(cleanup_intermediate_bcf)
    
    script:
    """
    log_msg() {
        local level="\$1"
        local msg="\$2"
        echo "[\$(date -Iseconds)] [${sample}] [CLEANUP_SAMPLE_TEMP] [\${level}] \${msg}"
    }

    log_msg "INFO" "Starting sample cleanup process (symlink-aware)"

    remove_file_follow_symlink() {
        local f="\$1"
        if [ -L "\${f}" ]; then
            local target
            target=\$(realpath "\${f}")
            if [ -e "\$target" ]; then
                rm -f "\$target" && log_msg "INFO" "Cleaned up symlink target: \$target"
            fi
            rm -f "\$f" && log_msg "INFO" "Cleaned up symlink: \$f"
        elif [ -e "\$f" ]; then
            rm -f "\$f" && log_msg "INFO" "Cleaned up file: \$f"
        fi
    }

    for vcf in ${output_vcf_snps}; do
        [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        [ -e "\${vcf}.csi" ] && remove_file_follow_symlink "\${vcf}.csi"
    done

    for vcf in ${output_vcf_indels}; do
        [ -e "\${vcf}" ] && remove_file_follow_symlink "\${vcf}"
        [ -e "\${vcf}.csi" ] && remove_file_follow_symlink "\${vcf}.csi"
    done

    if [ "${cleanup_intermediate_bam}" = "true" ]; then
        [ -e "${bam}" ] && remove_file_follow_symlink "${bam}"
        [ -e "${bam_index}" ] && remove_file_follow_symlink "${bam_index}"
        log_msg "INFO" "Cleaned up BAM files"
    fi

    [ -e "${coverage_bed}" ] && remove_file_follow_symlink "${coverage_bed}"
    [ -e "${zero_bcf}" ] && remove_file_follow_symlink "${zero_bcf}"

    if [ "${cleanup_intermediate_bcf}" = "true" ]; then
        for vcf in "${sample}.snps_"*.bcf "${sample}.indels_"*.bcf; do
            [ -e "\$vcf" ] && remove_file_follow_symlink "\$vcf"
        done
        log_msg "INFO" "Cleaned up consensus VCFs"
    fi

    for symlink in ./*; do
        if [ -L "\$symlink" ]; then
            rm -f "\$symlink" && log_msg "INFO" "Cleaned up symlink in current directory: \$symlink"
        fi
    done

    log_msg "INFO" "Sample cleanup completed"
    """
}
