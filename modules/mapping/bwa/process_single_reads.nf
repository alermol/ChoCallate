process MAP_BWA_SINGLE {
    maxForks 1
    cpus params.mapping.cpu
    //afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} "\$(realpath tmp/ref_genome.fasta)" tmp/read1 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """
}