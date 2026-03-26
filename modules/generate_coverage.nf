process GENERATE_COVERAGE {
    cpus 1
    maxForks 1
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")

    output:
    path("coverage.bed"), emit: coverage

    script:
    """
    samtools depth -J --threads ${task.cpus} tmp/input.bam | awk '\$3 >= ${params.coverage.min_coverage} {print \$1,\$2-1,\$2}' | bedops --merge - > coverage.bed
    """

    stub:
    """
    touch coverage.bed
    """
}

process INTERSECT_CUSTOM_BED {
    cpus 1
    maxForks params.coverage.forks
    afterScript 'stage_cleanup.sh'

    input:
    path("tmp/coverage.bed")
    path("tmp/custom.bed")

    output:
    path("coverage.bed"), emit: coverage

    script:
    """
    bedops --intersect tmp/coverage.bed tmp/custom.bed > coverage.bed
    """

    stub:
    """
    touch coverage.bed
    """
}