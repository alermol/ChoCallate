process MAP_BWA_PAIRED {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.bwa.cpu
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2")
    val(ref_index)

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    def rg_id = "@RG\\tID:${sample_id}\\tSM:${sample_id}"
    """
    bwa mem -t ${task.cpus} ${params.bam_preparation.mapping.bwa.extra_args} -R "${rg_id}" "${ref_index}" tmp/read1 tmp/read2 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}

process MAP_BWA_SINGLE {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.bwa.cpu
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1")
    val(ref_index)

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${params.bam_preparation.mapping.bwa.extra_args} "${ref_index}" tmp/read1 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}