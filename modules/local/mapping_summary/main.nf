process MAPPING_SUMMARY {
    label 'process_low'
    cpus 1
    tag "mapping_summary"

    publishDir(
        "${params.output.directory}/mapping_stats",
        mode: params.output.publish_mode,
        pattern: 'mapping_summary.tsv'
    )

    input:
    val sample_ids
    path(raw_flagstats, stageAs: 'raw/*')
    path(final_flagstats, stageAs: 'flagstat/*')
    path(stats_files, stageAs: 'stats/*')
    path(idxstats_files, stageAs: 'idxstats/*')

    output:
    path 'mapping_summary.tsv'

    script:
    """
    printf '%s\\n' ${sample_ids.join(' ')} > sample_order.txt

    mapping_stats.py \
        --sample-order sample_order.txt \
        --raw-dir raw \
        --flagstat-dir flagstat \
        --stats-dir stats \
        --idxstats-dir idxstats \
        --output mapping_summary.tsv
    """
}
