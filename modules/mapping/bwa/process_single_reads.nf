process MAP_BWA_SINGLE {
    maxForks 1
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    GENOME_BASENAME=\$(basename \$(realpath tmp/ref_genome.fasta))
    
    bwa mem ${params.mapping.extra_args} -t ${task.cpus} ${params.input.reference_index_dir}/\${GENOME_BASENAME} tmp/read1 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """
}