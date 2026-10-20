process FLAKY {
    label 'process_low'
    publishDir params.outdir, mode: 'copy'

    output:
    path 'ok.txt'
    path 'memory.txt'
    path 'time.txt'

    script:
    """
    if [ ${task.attempt} -lt 2 ]; then
        kill -9 \$\$
    fi
    printf '%s\\n' "${task.memory.toMega()}" > memory.txt
    printf '%s\\n' "${task.time.toMinutes()}" > time.txt
    touch ok.txt
    """
}

workflow {
    FLAKY()
}
