process CALLING_BCFTOOLS {
    label 'process_medium'

    script:
    """
    printf '%s\\n' "${task.attempt}" >> "${params.attempt_log}"
    exit ${params.caller_exit_status ?: 1}
    """
}

workflow {
    CALLING_BCFTOOLS()
}
