process MAP_MINIMAP2_PAIRED {
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
    GENOME_BASENAME=\$(basename \$(realpath tmp/ref_genome.fasta))
    
    minimap2 ${params.mapping.extra_args} -ax sr -t ${task.cpus} -a ${params.input.reference_index_dir}/\${GENOME_BASENAME}.mmi tmp/read1 tmp/read2 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """
}