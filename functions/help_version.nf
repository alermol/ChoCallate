def getPipelineMetadata() {
    return [
        name: 'ChoCallate',
        version: '1.0.0',
        author: 'A. Ermolaev',
        description: 'Consensus variant calling pipeline',
        license: 'MIT',
        homepage: 'https://github.com/alermol/chocallate'
    ]
}

def showVersion() {
    def meta = getPipelineMetadata()
    
    println "${meta.name} - ${meta.description}"
    println "Version: ${meta.version}"
    println "Author:  ${meta.author}"
    println "License: ${meta.license}"
    println "Home:    ${meta.homepage}"
}

def showHelp() {
    def meta = getPipelineMetadata()
    
    println "${meta.name} - ${meta.description}"
    println "Version: ${meta.version}\n"
    println ""
    println "DESCRIPTION:"
    println "${meta.description}\n"
    println ""
    println "USAGE:"
    println "  nextflow run main.nf [OPTIONS]\n"
    println ""
    println "BASIC OPTIONS:"
    println "  --reference_genome <file>    Reference genome FASTA file"
    println "  --reference_index <dir>      Bowtie2 index directory (for FASTQ input)"
    println "  --samples_tsv <file>         Samples TSV file"
    println "  --outdir <dir>               Output directory\n"
    println ""
    println "INPUT FORMAT OPTIONS:"
    println "  --input_format <format>      Input format: 'fastq' or 'bam' (default: fastq)"
    println "  --reads_type <type>          Read type: 'se', 'pe', or 'mx' (default: pe)"
    println "  --reads_source <source>      Source: 'gbs' or 'wgs' (default: gbs)\n"
    println ""
    println "VARIANT CALLING OPTIONS:"
    println "  --ploidy <number>            Ploidy level (default: 2)"
    println "  --effective_callers <list>   Comma-separated caller list or '-' for auto"
    println "                               Available: bcftools,gatk,freebayes,snver,vardict\n"
    println ""
    println "QUALITY FILTERING:"
    println "  --min_coverage <number>      Minimum coverage (default: 5)"
    println "  --min_base_quality <number>  Minimum base quality (default: 5)"
    println "  --min_map_qual <number>      Minimum mapping quality (default: 5)"
    println "  --min_snp_qual <number>      Minimum SNP quality (default: 5)\n"
    println ""
    println "CONSENSUS OPTIONS:"
    println "  --cons_type <type>           Consensus type: 'mj', 'n1', or 'fc' (default: mj)"
    println "                               mj=majority, n1=n-1, fc=full consensus"
    println "  --win_size <number>          Window size for consensus (default: 1000000)\n"
    println ""
    println "RESOURCE ALLOCATION:"
    println "  --bowtie2_cpu <number>       CPUs for Bowtie2 (default: 10)"
    println "  --cons_cpus <number>         CPUs for consensus (default: 5)\n"
    println ""
    println "OTHER OPTIONS:"
    println "  --debug                      Enable debug mode (preserves intermediate files)"
    println "  --version                    Show version information"
    println "  --help                       Show this help message\n"
    println ""
    println "EXAMPLES:"
    println "  # Basic FASTQ input:"
    println "  nextflow run main.nf --reference_genome ref.fasta --reference_index ref_idx \\"
    println "                       --samples_tsv samples.tsv --outdir results\n"
    println ""
    println "  # BAM input:"
    println "  nextflow run main.nf --reference_genome ref.fasta --input_format bam \\"
    println "                       --samples_tsv samples_bam.tsv --outdir results\n"
    println ""
    println "  # Polyploid with specific callers:"
    println "  nextflow run main.nf --reference_genome ref.fasta --reference_index ref_idx \\"
    println "                       --samples_tsv samples.tsv --ploidy 4 \\"
    println "                       --effective_callers gatk,freebayes,snver\n"
    println ""
    println "For more information, visit: ${meta.homepage}"
}

def handleHelpVersion() {
    if (params.version) {
        showVersion()
        System.exit(0)
    }
    
    if (params.help) {
        showHelp()
        System.exit(0)
    }
}
