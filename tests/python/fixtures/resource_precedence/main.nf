process MAPPING_BOWTIE2 {
    label 'process_medium'
    cpus 99
    publishDir params.outdir, mode: 'copy'

    output:
    path 'resources.txt'

    script:
    """
    printf 'cpus=%s\\n' "${task.cpus}" > resources.txt
    printf 'memory=%s\\n' "${task.memory.toMega()}" >> resources.txt
    """
}

workflow {
    MAPPING_BOWTIE2()
}
