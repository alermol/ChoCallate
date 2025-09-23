def getAvailableCallersCount(String available_callers) {
    def callersList = available_callers.split(/\s*,\s*/).findAll { it }
    return callersList.size()
}

def getConsensusThreshold(String cons_type, String available_callers) {
    if (cons_type == 'mj') {
        return getAvailableCallersCount(available_callers).intdiv(2) + 1
    } else if (cons_type == 'n1') {
        return getAvailableCallersCount(available_callers) - 1
    } else if (cons_type == 'fc') {
        return getAvailableCallersCount(available_callers)
    } else {
        error "Invalid consensus type: ${cons_type}. Available types: mj, n1, fc"
    }
}

def normalizeCallerNames(String callers) {
    def callerList = callers.split(/\s*,\s*/).findAll { it }
    return callerList.collect { it.toLowerCase().trim() }.join(',')
}

def allEffectiveCallersInAvailable(String effective_callers, String available_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }.collect { it.toLowerCase().trim() }
    def availList = available_callers.split(/\s*,\s*/).findAll { it }.collect { it.toLowerCase().trim() }
    def missing = effList.findAll { !(it in availList) }
    if (missing) {
        println "The following effective callers are not in available callers: ${missing.join(', ')}"
        return false
    }
    return true
}

def effectiveCallersAtLeastThree(String effective_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }
    return effList.size()
}

def allEffectiveCallersDiploidSuitable(String effective_callers, String diploid_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }.collect { it.toLowerCase().trim() }
    def dipList = diploid_callers.split(/\s*,\s*/).findAll { it }.collect { it.toLowerCase().trim() }
    def unsuitable = effList.findAll { !(it in dipList) }
    if (unsuitable) {
        println "The following effective callers are not suitable for diploid calling: ${unsuitable.join(', ')}"
        return false
    }
    return true
}

def allEffectiveCallersPolyploidSuitable(String effective_callers, String polyploid_callers) {
    def effList = effective_callers.split(/\s*,\s*/).findAll { it }.collect { it.toLowerCase().trim() }
    def polyList = polyploid_callers.split(/\s*,\s*/).findAll { it }.collect { it.toLowerCase().trim() }
    def unsuitable = effList.findAll { !(it in polyList) }
    if (unsuitable) {
        println "The following effective callers are not suitable for polyploid calling: ${unsuitable.join(', ')}"
        return false
    }
    return true
}


def validateFile(String filePath, String fileType) {
    if (!filePath) {
        return [valid: false, error: "${fileType} path is not specified"]
    }
    
    def file = new File(filePath)
    if (!file.exists()) {
        return [valid: false, error: "${fileType} file does not exist: ${filePath}"]
    }
    
    if (!file.canRead()) {
        return [valid: false, error: "${fileType} file is not readable: ${filePath}"]
    }
    
    if (file.length() == 0) {
        return [valid: false, error: "${fileType} file is empty: ${filePath}"]
    }
    
    return [valid: true, file: file]
}

def validateInputFormat(String inputFormat) {
    def valid = ['fastq','bam']
    if (!inputFormat || !valid.contains(inputFormat.toLowerCase())) {
        return [valid: false, error: "Invalid input_format: ${inputFormat}. Valid options: fastq, bam"]
    }
    return [valid: true, value: inputFormat.toLowerCase()]
}

def validateTSVFile(String filePath, String inputFormat, String readsType) {
    def fileValidation = validateFile(filePath, "Samples TSV")
    if (!fileValidation.valid) {
        return fileValidation
    }
    
    def file = fileValidation.file
    // Read file content as string first to preserve tab characters
    def content = file.text
    def lines = content.split('\n')
    
    if (lines.size() < 1) {
        return [valid: false, error: "Samples TSV file must contain at least 1 sample line"]
    }
    
    // No header - all lines are sample data
    def startIndex = 0
    
    // Normalize for comparisons
    def fmt = (inputFormat ?: 'fastq').toLowerCase()
    def rt = (readsType ?: '').toLowerCase()

    // Helpers to support comma-separated lists of paths in cells
    def parsePaths = { String s ->
        (s ?: '')
            .split(/\s*,\s*/)
            .collect { it.trim() }
            .findAll { it }
    }

    def anyMissingPath = { List<String> paths ->
        for (p in paths) {
            if (!new File(p).exists()) {
                return p
            }
        }
        return null
    }

    // Validate sample lines
    for (int i = startIndex; i < lines.size(); i++) {
        def line = lines[i].trim()
        if (line.isEmpty()) continue
        
        def columns = line.split('\t')
        if (fmt == 'fastq') {
            if (columns.size() < 4) {
                return [valid: false, error: "Sample line ${i+1} has insufficient columns: ${columns.size()} < 4. Expected 4 columns for FASTQ mode"]
            }
            if (!columns[0] || columns[0].trim().isEmpty()) {
                return [valid: false, error: "Sample line ${i+1} has empty sample ID"]
            }
            if (rt == 'se') {
                def sePaths = parsePaths(columns[3])
                if (sePaths.isEmpty()) {
                    return [valid: false, error: "FASTQ se mode requires column 4 (SE read) in line ${i+1}"]
                }
                def missing = anyMissingPath(sePaths)
                if (missing) {
                    return [valid: false, error: "SE FASTQ file listed in column 4 does not exist (line ${i+1}): ${missing}"]
                }
            } else if (rt == 'pe') {
                def r1Paths = parsePaths(columns[1])
                def r2Paths = parsePaths(columns[2])
                if (r1Paths.isEmpty() || r2Paths.isEmpty()) {
                    return [valid: false, error: "FASTQ pe mode requires columns 2 and 3 (R1,R2) in line ${i+1}"]
                }
                if (r1Paths.size() != r2Paths.size()) {
                    return [valid: false, error: "FASTQ pe mode requires equal number of R1 and R2 files (line ${i+1}): R1=${r1Paths.size()}, R2=${r2Paths.size()}"]
                }
                def missingR1 = anyMissingPath(r1Paths)
                if (missingR1) {
                    return [valid: false, error: "PE FASTQ R1 listed in column 2 does not exist (line ${i+1}): ${missingR1}"]
                }
                def missingR2 = anyMissingPath(r2Paths)
                if (missingR2) {
                    return [valid: false, error: "PE FASTQ R2 listed in column 3 does not exist (line ${i+1}): ${missingR2}"]
                }
            } else if (rt == 'mx') {
                def r1Paths = parsePaths(columns[1])
                def r2Paths = parsePaths(columns[2])
                def sePaths = parsePaths(columns[3])
                if (r1Paths.isEmpty() || r2Paths.isEmpty() || sePaths.isEmpty()) {
                    return [valid: false, error: "FASTQ mx mode requires columns 2,3,4 (R1,R2,SE) in line ${i+1}"]
                }
                if (r1Paths.size() != r2Paths.size()) {
                    return [valid: false, error: "FASTQ mx mode requires equal number of R1 and R2 files (line ${i+1}): R1=${r1Paths.size()}, R2=${r2Paths.size()}"]
                }
                def missingR1 = anyMissingPath(r1Paths)
                if (missingR1) return [valid: false, error: "MX FASTQ R1 listed in column 2 does not exist (line ${i+1}): ${missingR1}"]
                def missingR2 = anyMissingPath(r2Paths)
                if (missingR2) return [valid: false, error: "MX FASTQ R2 listed in column 3 does not exist (line ${i+1}): ${missingR2}"]
                def missingSE = anyMissingPath(sePaths)
                if (missingSE) return [valid: false, error: "MX FASTQ SE listed in column 4 does not exist (line ${i+1}): ${missingSE}"]
            } else {
                return [valid: false, error: "Invalid reads_type '${readsType}' for FASTQ input. Expected one of: se, pe, mx"]
            }
        } else if (fmt == 'bam') {
            if (columns.size() < 2) {
                return [valid: false, error: "Sample line ${i+1} has insufficient columns: ${columns.size()} < 2. Expected: sample_id, bam_path"]
            }
            if (!columns[0] || columns[0].trim().isEmpty()) {
                return [valid: false, error: "Sample line ${i+1} has empty sample ID"]
            }
            def bamPath = columns[1]
            def bamPaths = parsePaths(bamPath)
            if (bamPaths.isEmpty()) {
                return [valid: false, error: "BAM path is required in column 2 (line ${i+1})"]
            }
            for (p in bamPaths) {
                def bf = new File(p)
                if (!bf.exists()) {
                    return [valid: false, error: "BAM file listed in column 2 does not exist (line ${i+1}): ${p}"]
                }
                if (!p.toLowerCase().endsWith('.bam')) {
                    return [valid: false, error: "BAM file in line ${i+1} must have .bam extension: ${p}"]
                }
            }
        }
    }
    
    def sampleCount = lines.size() - startIndex
    if (sampleCount < 1) {
        return [valid: false, error: "No valid sample lines found in TSV file"]
    }
    
    return [valid: true, file: file, sampleCount: sampleCount]
}

def validateReferenceGenome(String filePath) {
    def fileValidation = validateFile(filePath, "Reference genome")
    if (!fileValidation.valid) {
        return fileValidation
    }
    
    def file = fileValidation.file
    def fileName = file.getName().toLowerCase()
    
    // Check file extension (including gzipped versions)
    if (!fileName.endsWith('.fasta') && !fileName.endsWith('.fa') && !fileName.endsWith('.fna') && 
        !fileName.endsWith('.fasta.gz') && !fileName.endsWith('.fa.gz') && !fileName.endsWith('.fna.gz')) {
        return [valid: false, error: "Reference genome must be a FASTA file (.fasta, .fa, .fna) or gzipped FASTA (.fasta.gz, .fa.gz, .fna.gz): ${filePath}"]
    }
    
    // Check if it's a valid FASTA file
    def firstLine
    if (fileName.endsWith('.gz')) {
        // For gzipped files, use gzip input stream
        firstLine = file.withInputStream { inputStream ->
            new java.util.zip.GZIPInputStream(inputStream).withReader { gzReader ->
                gzReader.readLine()
            }
        }
    } else {
        // For uncompressed files, read directly
        firstLine = file.withReader { reader ->
            reader.readLine()
        }
    }
    
    if (!firstLine || !firstLine.startsWith('>')) {
        return [valid: false, error: "Reference genome file is not a valid FASTA file (does not start with '>'): ${filePath}"]
    }
    
    return [valid: true, file: file]
}

def validateBowtie2Index(String indexPath) {
	if (!indexPath) {
		return [valid: false, error: "Bowtie2 index path is not specified"]
	}

	def pathObj = new File(indexPath)
	def indexDir = pathObj.isDirectory() ? pathObj : pathObj.getParentFile()
	if (!indexDir || !indexDir.exists()) {
		return [valid: false, error: "Bowtie2 index directory does not exist: ${indexPath}"]
	}

	// Try to infer the index prefix from the provided path
	def name = pathObj.isDirectory() ? null : pathObj.getName()
	def prefixGuess = null
	if (name) {
		// Remove known suffix patterns to get the prefix base
		prefixGuess = name.replaceAll(/(\.rev\.(1|2)\.(bt2|bt2l)|\.(1|2|3|4)\.(bt2|bt2l))$/, '')
		if (!prefixGuess || prefixGuess == name) {
			// If no suffix matched, assume the provided name is already the prefix
			prefixGuess = name
		}
	}

	// Helper to check existence of forward index set for a given prefix and extension
	def forwardSetExists = { String prefix, String ext ->
		[1,2,3,4].every { n -> new File(indexDir, "${prefix}.${n}.${ext}").exists() }
	}

	def found = false
	def usedPrefix = null
	def usedExt = null

	if (prefixGuess) {
		// Check for either bt2l or bt2 forward sets
		for (ext in ['bt2l','bt2']) {
			if (forwardSetExists(prefixGuess, ext)) {
				found = true
				usedPrefix = prefixGuess
				usedExt = ext
				break
			}
		}
	}

	// If not found and no explicit prefix provided (directory or non-matching file),
	// attempt to discover any valid prefix in the directory by scanning files
	if (!found) {
		def files = indexDir.listFiles() ?: []
		// Collect candidate prefixes from files matching .1.bt2(l)
		def candidates = files.collect { f -> f.getName() }
			.findAll { it ==~ /.+\.1\.(bt2|bt2l)$/ }
			.collect { it.replaceAll(/\.1\.(bt2|bt2l)$/, '') }
		// Check candidates
		for (cand in candidates) {
			for (ext in ['bt2l','bt2']) {
				if (forwardSetExists(cand, ext)) {
					found = true
					usedPrefix = cand
					usedExt = ext
					break
				}
			}
			if (found) break
		}
	}

	if (!found) {
		return [
			valid: false,
			error: "Bowtie2 index validation failed in '${indexDir}'. Expected forward index files '.1..4' with extension '.bt2' or '.bt2l' for a common prefix. Reverse index files are optional."]
	}

	// Build the list of discovered forward index files for reporting
	def indexFiles = [1,2,3,4].collect { n -> new File(indexDir, "${usedPrefix}.${n}.${usedExt}") }
	return [valid: true, indexDir: indexDir, indexFiles: indexFiles, prefix: usedPrefix, ext: usedExt]
}

def validateNumericParameter(Number value, String paramName, Number minValue, Number maxValue = null) {
    if (value == null) {
        return [valid: false, error: "${paramName} is not specified"]
    }
    
    if (value < minValue) {
        return [valid: false, error: "${paramName} (${value}) is below minimum value (${minValue})"]
    }
    
    if (maxValue != null && value > maxValue) {
        return [valid: false, error: "${paramName} (${value}) is above maximum value (${maxValue})"]
    }
    
    return [valid: true, value: value]
}

def validateCPUParameter(Number value, String paramName) {
    def maxCPUs = Runtime.runtime.availableProcessors()
    return validateNumericParameter(value, paramName, 1, maxCPUs)
}

def validateForkParameter(Number value, String paramName) {
    def maxCPUs = Runtime.runtime.availableProcessors()
    def maxForks = (maxCPUs > 1) ? (maxCPUs - 1) : 1
    return validateNumericParameter(value, paramName, 1, maxForks)
}

def validateQualityParameter(Number value, String paramName, Number maxValue) {
    return validateNumericParameter(value, paramName, 0, maxValue) // Quality scores typically 0-maxValue
}

def validatePloidy(Number ploidy) {
    if (ploidy == null) {
        return [valid: false, error: "Ploidy is not specified"]
    }
    
    if (ploidy < 2) {
        return [valid: false, error: "Ploidy (${ploidy}) must be 2 or greater"]
    }
    
    return [valid: true, value: ploidy]
}

def validateReadsType(String readsType) {
    def validTypes = ['se', 'pe', 'mx']
    if (!readsType || !validTypes.contains(readsType)) {
        return [valid: false, error: "Invalid reads_type: ${readsType}. Valid options: ${validTypes.join(', ')}"]
    }
    return [valid: true, value: readsType]
}

def validateReadsSource(String readsSource) {
    def validSources = ['gbs', 'wgs']
    if (!readsSource || !validSources.contains(readsSource)) {
        return [valid: false, error: "Invalid reads_source: ${readsSource}. Valid options: ${validSources.join(', ')}"]
    }
    return [valid: true, value: readsSource]
}

def validateConsensusType(String consType) {
    def validTypes = ['mj', 'n1', 'fc']
    if (!consType || !validTypes.contains(consType)) {
        return [valid: false, error: "Invalid cons_type: ${consType}. Valid options: ${validTypes.join(', ')}"]
    }
    return [valid: true, value: consType]
}

def validateLogLevel(String logLevel) {
    def validLevels = ['DEBUG', 'INFO', 'WARN', 'ERROR', 'FATAL']
    if (!logLevel || !validLevels.contains(logLevel.toUpperCase())) {
        return [valid: false, error: "Invalid log_level: ${logLevel}. Valid options: ${validLevels.join(', ')}"]
    }
    return [valid: true, value: logLevel.toUpperCase()]
}

def validateLogFormat(String logFormat) {
    def validFormats = ['json', 'text', 'both']
    if (!logFormat || !validFormats.contains(logFormat.toLowerCase())) {
        return [valid: false, error: "Invalid log_format: ${logFormat}. Valid options: ${validFormats.join(', ')}"]
    }
    return [valid: true, value: logFormat.toLowerCase()]
}

def validateWindowSize(Number winSize) {
    if (winSize == null) {
        return [valid: false, error: "Window size is not specified"]
    }
    
    if (winSize < 1000000) {
        return [valid: false, error: "Window size (${winSize}) is too small. Minimum: 1000000 bp"]
    }
    
    if (winSize > 10000000) {
        return [valid: false, error: "Window size (${winSize}) is too large. Maximum: 10000000 bp"]
    }
    
    return [valid: true, value: winSize]
}



def validateAllParameters(Map params) {
    def errors = []
    def warnings = []
    
    // Input format
    def inputFormatValidation = validateInputFormat(params.input_format)
    if (!inputFormatValidation.valid) {
        errors << inputFormatValidation.error
    }
    def inputFormat = inputFormatValidation.valid ? inputFormatValidation.value : 'fastq'

    // File validations
    def samplesValidation = validateTSVFile(params.samples_tsv, inputFormat, params.reads_type)
    if (!samplesValidation.valid) {
        errors << samplesValidation.error
    }
    
    def refGenomeValidation = validateReferenceGenome(params.reference_genome)
    if (!refGenomeValidation.valid) {
        errors << refGenomeValidation.error
    }
    
    if (inputFormat == 'fastq') {
        def refIndexValidation = validateBowtie2Index(params.reference_index)
        if (!refIndexValidation.valid) {
            errors << refIndexValidation.error
        }
    }
    
    // Numeric parameter validations
    def ploidyValidation = validatePloidy(params.ploidy)
    if (!ploidyValidation.valid) {
        errors << ploidyValidation.error
    }
    
    def minCoverageValidation = validateQualityParameter(params.min_coverage, "min_coverage", 10000000)
    if (!minCoverageValidation.valid) {
        errors << minCoverageValidation.error
    }
    
    def minBaseQualityValidation = validateQualityParameter(params.min_base_quality, "min_base_quality", 10000000)
    if (!minBaseQualityValidation.valid) {
        errors << minBaseQualityValidation.error
    }
    
    def minMapQualValidation = validateQualityParameter(params.min_map_qual, "min_map_qual", 10000000)
    if (!minMapQualValidation.valid) {
        errors << minMapQualValidation.error
    }
    
    def minSnpQualValidation = validateQualityParameter(params.min_snp_qual, "min_snp_qual", 10000000)
    if (!minSnpQualValidation.valid) {
        errors << minSnpQualValidation.error
    }
    
    // CPU and resource validations
    if (inputFormat == 'fastq') {
        def bowtie2CpuValidation = validateCPUParameter(params.bowtie2_cpu, "bowtie2_cpu")
        if (!bowtie2CpuValidation.valid) {
            errors << bowtie2CpuValidation.error
        }
    }
    
    def bcftoolsCpuValidation = validateCPUParameter(params.bcftools_cpu, "bcftools_cpu")
    if (!bcftoolsCpuValidation.valid) {
        errors << bcftoolsCpuValidation.error
    }
    
    def vardictCpuValidation = validateCPUParameter(params.vardict_cpu, "vardict_cpu")
    if (!vardictCpuValidation.valid) {
        errors << vardictCpuValidation.error
    }
    
    def zeroBcfCpuValidation = validateCPUParameter(params.zero_bcf_cpu, "zero_bcf_cpu")
    if (!zeroBcfCpuValidation.valid) {
        errors << zeroBcfCpuValidation.error
    }
    
    def consCpusValidation = validateCPUParameter(params.cons_cpus, "cons_cpus")
    if (!consCpusValidation.valid) {
        errors << consCpusValidation.error
    }
    
    // Fork validations
    if (inputFormat == 'fastq') {
        def bowtie2ForksValidation = validateForkParameter(params.bowtie2_forks, "bowtie2_forks")
        if (!bowtie2ForksValidation.valid) {
            errors << bowtie2ForksValidation.error
        }
    }
    
    def callingForksValidation = validateForkParameter(params.calling_forks, "calling_forks")
    if (!callingForksValidation.valid) {
        errors << callingForksValidation.error
    }
    
    def zeroBcfForksValidation = validateForkParameter(params.zero_bcf_forks, "zero_bcf_forks")
    if (!zeroBcfForksValidation.valid) {
        errors << zeroBcfForksValidation.error
    }
    
    def consForksValidation = validateForkParameter(params.cons_forks, "cons_forks")
    if (!consForksValidation.valid) {
        errors << consForksValidation.error
    }
    
    // String parameter validations
    if (inputFormat == 'fastq') {
        def readsTypeValidation = validateReadsType(params.reads_type)
        if (!readsTypeValidation.valid) {
            errors << readsTypeValidation.error
        }
    } else {
        if (params.reads_type && params.reads_type !in ['se','pe','mx']) {
            warnings << "reads_type (${params.reads_type}) is ignored when input_format=bam"
        }
    }
    
    def readsSourceValidation = validateReadsSource(params.reads_source)
    if (!readsSourceValidation.valid) {
        errors << readsSourceValidation.error
    }
    
    def consTypeValidation = validateConsensusType(params.cons_type)
    if (!consTypeValidation.valid) {
        errors << consTypeValidation.error
    }
    
    // Other numeric validations
    def winSizeValidation = validateWindowSize(params.win_size)
    if (!winSizeValidation.valid) {
        errors << winSizeValidation.error
    }
    
    // Logging parameter validations
    if (params.log_level) {
        def logLevelValidation = validateLogLevel(params.log_level)
        if (!logLevelValidation.valid) {
            errors << logLevelValidation.error
        }
    }
    
    if (params.log_format) {
        def logFormatValidation = validateLogFormat(params.log_format)
        if (!logFormatValidation.valid) {
            errors << logFormatValidation.error
        }
    }
    
    // Warnings for potentially problematic configurations
    if (params.min_coverage < 2) {
        warnings << "min_coverage (${params.min_coverage}) is very low. This may result in many false positive variants."
    }
    
    if (params.min_base_quality < 5) {
        warnings << "min_base_quality (${params.min_base_quality}) is low. This may result in poor quality variant calls."
    }
    
    if (params.win_size > 5000000) {
        warnings << "win_size (${params.win_size}) is very large. This may impact memory usage and processing time."
    }
    
    return [valid: errors.isEmpty(), errors: errors, warnings: warnings]
}

