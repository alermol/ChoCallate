process PREPARE_BAM {
    maxForks params.bowtie2_forks
    cpus params.bowtie2_cpu

    tag "${sample_id}"

    input:
    tuple val(sample_id), val(read1), val(read2), val(read3)
    path(genome_dictionary)
    tuple path(ref_genome), path(genome_fai)
    val(ref_index)
    val(bowtie2_extra_args)

    output:
    tuple path("${sample_id}.bam"), path("${sample_id}.bam.csi"), emit: bam

    script:
    """
    prepare_bam.sh \
        ${params.input_format} \
        ${params.reads_type} \
        ${sample_id} \
        "${read1}" \
        "${read2}" \
        "${read3}" \
        "${ref_genome}" \
        "${genome_dictionary}" \
        "${genome_fai}" \
        "${ref_index}" \
        ${params.min_map_qual} \
        ${task.cpus} \
        "${params.cleanup_intermediate_subfolders}" \
        "${params.cleanup_input_symlinks}" \
        "${params.bowtie2_extra_args}"
    """
}


