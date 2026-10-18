process MAPPING_MINIMAP2 {
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
    // The aligner and BAM post-processing run concurrently in one pipe, so
    // split the requested CPUs. samtools view between them adds one thread.
    def cpus = task.cpus as int
    def aligner_cpus = Math.max(1, (cpus + 1).intdiv(2))
    def postprocess_cpus = Math.max(1, cpus - aligner_cpus)
    def memory_mb = task.memory ? task.memory.toMega() : 0
    def memory_option = memory_mb > 0 ? "--memory-mb ${memory_mb}" : ''
    def reads = params.input.reads_type == 'se' ? "tmp/read1" :
                params.input.reads_type == 'pe' ? "tmp/read1 tmp/read2" : ''
    """
    set -o pipefail

    minimap2 \
        ${params.mapping.extra_args} \
        -ax sr \
        -t ${aligner_cpus} \
        -R "${rg_id}" \
        -a "index/${index_prefix}.mmi" \
        ${reads} | \
    samtools view -u | \
    postprocess_bam.sh \
        --input - \
        --output mapping.bam \
        --reference tmp/ref_genome.fasta \
        --dict tmp/ref_genome.dict \
        --mapq ${params.bam_filter.min_map_qual} \
        --threads ${postprocess_cpus} \
        ${memory_option} \
        --java-options "${java_options}" \
        ${rmdup} \
        ${left_align}
    """
}
