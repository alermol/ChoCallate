process MAP_BWA_PAIRED {
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
    def rg_id = "@RG\\tID:${sample_id}\\tSM:${sample_id}"
    """
    GENOME_BASENAME=\$(basename \$(realpath tmp/ref_genome.fasta))
    
    bwa mem ${params.mapping.extra_args} -t ${task.cpus} -R "${rg_id}" ${params.input.reference_index_dir}/\${GENOME_BASENAME} tmp/read1 tmp/read2 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """
}