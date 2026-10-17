process MAPPING_BOWTIE2 {
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
    path("tmp/ref_genome_real.fasta")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.dict")
    path("tmp/ref_genome.fasta.fai")

    output:
    tuple val(sample_id), path("mapping.bam"),
        path("mapping.bam.csi"), emit: bam

    script:
    def rmdup = params.rmdup.enabled ? '--rmdup' : ''
    def left_align = params.left_align_indels.enabled ? '--left-align' : ''
    def platform = params.mapping.platform ?: 'ILLUMINA'
    def library = params.mapping.library ?: sample_id
    def index_prefix = params.input.reference_index_prefix
    def reads = params.input.reads_type == 'se' ? "-U tmp/read1" :
                params.input.reads_type == 'pe' ? "-1 tmp/read1 -2 tmp/read2" :
                params.input.reads_type == 'mx' ? "-1 tmp/read1 -2 tmp/read2 -U tmp/read3" : ''
    """
    set -o pipefail

    if [ -n "${index_prefix ?: ''}" ]; then
        INDEX_PREFIX="${index_prefix}"
    else
        INDEX_PREFIX="\$(basename \$(realpath tmp/ref_genome_real.fasta))"
    fi

    bowtie2 \
        ${params.mapping.extra_args} \
        --threads ${task.cpus} \
        --rg-id "${sample_id}" \
        --rg "SM:${sample_id}" \
        --rg "PL:${platform}" \
        --rg "LB:${library}" \
        -x "${params.input.reference_index_dir}/\$INDEX_PREFIX" \
        ${reads} | \
    samtools view --threads ${task.cpus} -b | \
    postprocess_bam.sh \
        --input - \
        --output mapping.bam \
        --reference tmp/ref_genome.fasta \
        --dict tmp/ref_genome.dict \
        --mapq ${params.bam_filter.min_map_qual} \
        --threads ${task.cpus} \
        ${rmdup} \
        ${left_align}
    """
}
