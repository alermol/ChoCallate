process GENERATE_COVERAGE {
    cpus 1
    maxForks 1
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path(include_bed)
    path(exclude_bed)

    output:
    path("coverage.bed"), emit: coverage

    script:
    def ibed = include_bed.name != "NO_FILE" ? "-b ${include_bed}" : ""
    def ebed = exclude_bed.name != "NO_FILE" ? "| bedops --difference - ${exclude_bed}" : ""
    """
    samtools depth -J ${ibed} --threads ${task.cpus} tmp/input.bam | awk '\$3 >= ${params.coverage.min_coverage} {print \$1,\$2-1,\$2}' ${ebed} | bedops --merge - > coverage.bed
    """
}
