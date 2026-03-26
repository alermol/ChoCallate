process FILTER_MAPPING_BAM {
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
    samtools view --threads "${task.cpus}" -b -q "${params.bam_filter.min_map_qual}" tmp/input.bam | samtools sort -n - | samtools fixmate --threads "${task.cpus}" -m - - | samtools sort --threads "${task.cpus}" - | samtools addreplacerg -w -r "@RG\tID:${sample_id}\tSM:${sample_id}" --threads "${task.cpus}" - -o output.bam
    """

    stub:
    """
    touch output.bam
    """
}

process FILTER_INPUT_BAM {
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
    samtools view -F 4 --threads "${task.cpus}" -b -q "${params.bam_filter.min_map_qual}" tmp/input.bam | samtools sort -n - | samtools fixmate --threads "${task.cpus}" -m - - | samtools sort --threads "${task.cpus}" - | samtools addreplacerg -w -r "@RG\tID:${sample_id}\tSM:${sample_id}" --threads "${task.cpus}" - -o output.bam
    """

    stub:
    """
    touch output.bam
    """
}