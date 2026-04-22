process MAP_BOWTIE2_MIXED {
    maxForks 1
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    """
    GENOME_BASENAME=\$(basename \$(realpath tmp/ref_genome.fasta))
    bowtie2 ${params.mapping.extra_args} --threads ${task.cpus} --rg-id "${sample_id}" --rg "SM:${sample_id}" -x ${params.input.reference_index_dir}/\${GENOME_BASENAME} -1 tmp/read1 -2 tmp/read2 -U tmp/read3 | samtools view --threads ${task.cpus} -b -o mapping.bam
    """
}
