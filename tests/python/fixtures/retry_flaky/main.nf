process FLAKY {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    output:
    path 'ok.txt'
    path 'memory.txt'

    script:
    """
    if [ ${task.attempt} -lt 2 ]; then
        kill -9 \$\$
    fi
    printf '%s\\n' "${task.memory.toMega()}" > memory.txt
    touch ok.txt
    """
}

workflow {
    FLAKY()
}
