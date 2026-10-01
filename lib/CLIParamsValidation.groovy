class CLIParamsValidation {
    static void reference_genome_validation(String reference_genome) {
        if (reference_genome == null) {
            println "ERROR: Reference genome is required"
            System.exit(1)
        }
    }

    static void samples_tsv_validation(String samples_tsv) {
        if (samples_tsv == null) {
            println "ERROR: Samples TSV is required"
            System.exit(1)
        }
    }

    static void reference_index_dir_validation(String reference_index_dir) {
        if (reference_index_dir == null) {
            println "ERROR: Reference index directory is required"
            System.exit(1)
        }
    }

    static void effective_callers_validation(String callers) {
        if (callers == null) {
            println "ERROR: Effective callers are not specified"
            System.exit(1)
        }
        def available_callers = ['bcftools', 'gatk', 'freebayes']
        def callersList = callers.split(',').collect { it.trim() }
        def unknown_callers = callersList.findAll { !(it in available_callers) }
        if (unknown_callers) {
            println "ERROR: Unknown caller(s): ${unknown_callers}. Available callers are: ${available_callers}"
            System.exit(1)
        }
    }

    static void mapper_validation(String mapper, String reads_type) {
        def available_mappers = ['bowtie2', 'bwa', 'minimap2']
        if (mapper == null) {
            println "ERROR: Mapper is required"
            System.exit(1)
        }
        if (!(mapper in available_mappers)) {
            println "ERROR: Invalid mapper: ${mapper}. Available mappers are: ${available_mappers}"
            System.exit(1)
        }
        if (mapper != 'bowtie2' && reads_type == 'mx') {
            println "ERROR: Mixed reads type is only supported for bowtie2 mapper"
            System.exit(1)
        }
    }

    static void cons_threshold_validation(Number cons_threshold, String callers) {
        if (cons_threshold == null) {
            println "ERROR: Consensus threshold is required"
            System.exit(1)
        }
        if (callers == null) {
            println "ERROR: Effective callers are not specified"
            System.exit(1)
        }
        def callersList = callers.split(',').collect { it.trim() }
        if (cons_threshold > callersList.size()) {
            println "ERROR: Consensus threshold must be less or equal to the number of effective callers"
            System.exit(1)
        }
        if (cons_threshold < 1) {
            println "ERROR: Consensus threshold must be greater than 0"
            System.exit(1)
        }
    }
}