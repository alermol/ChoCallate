process MAP_BOWTIE2_MIXED {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.bowtie2.cpu
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    bowtie2 --threads ${task.cpus} ${params.bam_preparation.mapping.bowtie2.extra_args} --rg-id "${sample_id}" --rg "SM:${sample_id}" -x "\$(realpath tmp/ref_genome.fasta)" -1 tmp/read1 -2 tmp/read2 -U tmp/read3 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}

process MAP_BOWTIE2_PAIRED {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.bowtie2.cpu
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    bowtie2 --threads ${task.cpus} ${params.bam_preparation.mapping.bowtie2.extra_args} -x "\$(realpath tmp/ref_genome.fasta)" -1 tmp/read1 -2 tmp/read2 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}

process MAP_BOWTIE2_SINGLE {
    maxForks params.bam_preparation.forks
    cpus params.bam_preparation.mapping.bowtie2.cpu
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    bowtie2 --threads ${task.cpus} ${params.bam_preparation.mapping.bowtie2.extra_args} -x "\$(realpath tmp/ref_genome.fasta)" -U tmp/read1 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """

    stub:
    """
    touch mapping.bam
    """
}