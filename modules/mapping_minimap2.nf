process MAP_MINIMAP2_PAIRED {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.minimap2.cpu
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2")
    val(ref_index)

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    minimap2 -ax sr -t ${task.cpus} ${params.bam_preparation.mapping.minimap2.extra_args} -a "${ref_index}" tmp/read1 tmp/read2 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}

process MAP_MINIMAP2_SINGLE {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.minimap2.cpu
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1")
    val(ref_index)

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    minimap2 -ax sr -t ${task.cpus} ${params.bam_preparation.mapping.minimap2.extra_args} "${ref_index}" tmp/read1 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}