process MAPPING_BWA {
    label 'process_medium'
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
    path(index_files, stageAs: 'index/*')
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
    def rg_id = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:${platform}\\tLB:${library}"
    def index_prefix = params.input.reference_index_prefix ?:
        file(params.input.reference_genome).name
    def java_options = "-Djava.io.tmpdir=\$TMPDIR " +
        ResourceUtils.javaHeapOption(task.memory)
    def reads = params.input.reads_type == 'se' ? "tmp/read1" :
                params.input.reads_type == 'pe' ? "tmp/read1 tmp/read2" : ''
    """
    set -o pipefail

    bwa mem \
        ${params.mapping.extra_args} \
        -t ${task.cpus} \
        -R "${rg_id}" \
        "index/${index_prefix}" \
        ${reads} | \
    samtools view --threads ${task.cpus} -b | \
    postprocess_bam.sh \
        --input - \
        --output mapping.bam \
        --reference tmp/ref_genome.fasta \
        --dict tmp/ref_genome.dict \
        --mapq ${params.bam_filter.min_map_qual} \
        --threads ${task.cpus} \
        --java-options "${java_options}" \
        ${rmdup} \
        ${left_align}
    """
}
