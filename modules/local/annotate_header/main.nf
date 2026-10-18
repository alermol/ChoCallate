process ANNOTATE_HEADER {
    // The header text is generated from NXF_CLI at task runtime, so a change
    // to the launch command line does not alter this process script.
    cache 'lenient'
    label 'process_low'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bcf")

    output:
    tuple val(sample_id), path("${sample_id}.bcf"), emit: consensus_bcf

    script:
    """
    write_command_header.py tmp/command_line.hdr
    bcftools annotate --threads ${task.cpus} \
        -h tmp/command_line.hdr \
        -Ob -o ${sample_id}.bcf tmp/input.bcf
    """
}
