process MAP_BOWTIE2_PAIRED {
    maxForks 1
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.dict")
    path("tmp/ref_genome.fasta.fai")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    def rmdup = params.rmdup.enabled ? '| picard MarkDuplicates -I - -O - -M /dev/null --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true' : ''
    def lai = params.left_align_indels.enabled ? '| gatk LeftAlignIndels -I /dev/stdin -O /dev/stdout -R tmp/ref_genome.fasta --sequence-dictionary tmp/ref_genome.dict -OBI false' : ''
    """
    GENOME_BASENAME=\$(basename \$(realpath tmp/ref_genome.fasta))
    
    bowtie2 \
        ${params.mapping.extra_args} \
        --threads ${task.cpus} \
        --rg-id ${sample_id} \
        --rg SM:${sample_id} \
        -x ${params.input.reference_index_dir}/\${GENOME_BASENAME} \
        -1 tmp/read1 \
        -2 tmp/read2 | \
    samtools view -F 4 --threads ${task.cpus} -b -q ${params.bam_filter.min_map_qual} | \
    samtools sort -n - | \
    samtools fixmate --threads ${task.cpus} -m - - | \
    samtools sort --threads ${task.cpus} ${rmdup} ${lai} | \
    samtools view -b -o mapping.bam
    """
}