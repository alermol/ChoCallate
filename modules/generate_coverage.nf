process GENERATE_COVERAGE {
    cpus 1
    maxForks 1
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path(custom_bed)

    output:
    path("coverage.bed"), emit: coverage

    script:
    def intersect_bed = custom_bed.name != "NO_FILE" ? "-b ${custom_bed}" : ""
    """
    samtools depth -J ${intersect_bed} --threads ${task.cpus} tmp/input.bam | awk '\$3 >= ${params.coverage.min_coverage} {print \$1,\$2-1,\$2}' | bedops --merge - > coverage.bed
    """
}
