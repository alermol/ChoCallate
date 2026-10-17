process PREPARE_BAM {
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
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
    """
    set -o pipefail

    samtools addreplacerg \
        -w \
        -r "ID:${sample_id}" \
        -r "SM:${sample_id}" \
        -r "PL:${platform}" \
        -r "LB:${library}" \
        -m overwrite_all \
        -u \
        -o - \
        tmp/read1 | \
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
