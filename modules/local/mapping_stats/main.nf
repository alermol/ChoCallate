process MAPPING_STATS {
    label 'process_low'
    cpus 1
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    // Raw aligner statistics, before MAPQ filtering and duplicate handling.
    publishDir(
        "${params.output.directory}/mapping_stats",
        mode: params.output.publish_mode,
        pattern: '*.raw.flagstat'
    )
    // Statistics of the BAM consumed by coverage and calling.
    publishDir(
        "${params.output.directory}/mapping_stats",
        mode: params.output.publish_mode,
        pattern: '*.final.flagstat'
    )
    publishDir(
        "${params.output.directory}/mapping_stats",
        mode: params.output.publish_mode,
        pattern: '*.stats'
    )
    publishDir(
        "${params.output.directory}/mapping_stats",
        mode: params.output.publish_mode,
        pattern: '*.idxstats'
    )

    input:
    tuple val(sample_id), path("tmp/input.bam"), path("tmp/input.bam.csi"),
        path("tmp/input.raw.flagstat")

    output:
    tuple val(sample_id),
        path("${sample_id}.raw.flagstat"),
        path("${sample_id}.final.flagstat"),
        path("${sample_id}.stats"),
        path("${sample_id}.idxstats"), emit: stats

    script:
    """
    cp tmp/input.raw.flagstat ${sample_id}.raw.flagstat

    samtools flagstat --threads ${task.cpus} \
        tmp/input.bam > ${sample_id}.final.flagstat
    samtools stats --threads ${task.cpus} \
        tmp/input.bam > ${sample_id}.stats
    samtools idxstats tmp/input.bam > ${sample_id}.idxstats
    """
}
