process MERGE_BEDS {
    cpus 1
    maxForks 1
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    input:
    path('tmp/?.bed', arity: '1..*')

    output:
    path('merged_coverage.bed'), emit: merged_coverage

    script:
    """
    bedops --merge tmp/*.bed > merged_coverage.bed
    """
}