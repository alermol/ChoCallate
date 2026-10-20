process REPORT_DROPPED_SAMPLES {
    label 'process_low'
    tag "dropped_samples"

    publishDir(
        params.output.directory,
        mode: params.output.publish_mode,
        pattern: 'dropped_samples.txt'
    )

    input:
    val dropped_samples

    output:
    path 'dropped_samples.txt'

    script:
    def lines = dropped_samples ? dropped_samples.join('\n') + '\n' : ''
    """
    printf '%s' '${lines}' > dropped_samples.txt
    """
}
