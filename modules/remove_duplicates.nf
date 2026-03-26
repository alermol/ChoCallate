process REMOVE_DUPLICATES {
    maxForks 1
    cpus 1
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")

    output:
    tuple val(sample_id), path("output.bam"), emit: bam

    script:
    """
    picard MarkDuplicates -I tmp/input.bam -O output.bam -M /dev/null --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true
    """

    stub:
    """
    touch output.bam
    """
}