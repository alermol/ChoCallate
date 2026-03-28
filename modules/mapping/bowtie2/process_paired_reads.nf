process MAP_BOWTIE2_PAIRED {
    maxForks 1
    cpus params.mapping.cpu
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    bowtie2 --threads ${task.cpus} -x "\$(realpath tmp/ref_genome.fasta)" -1 tmp/read1 -2 tmp/read2 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """
}