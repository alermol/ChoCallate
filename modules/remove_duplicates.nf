process REMOVE_DUPLICATES {
    maxForks params.bam_preparation.rm_duplicates.forks
    cpus params.bam_preparation.rm_duplicates.cpu
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")

    output:
    tuple val(sample_id), path("output.bam"), emit: bam

    script:
    """
    picard MarkDuplicates -I tmp/input.bam -O output.bam -M /dev/null ${params.bam_preparation.rm_duplicates.extra_args}
    """

    stub:
    """
    touch output.bam
    """
}