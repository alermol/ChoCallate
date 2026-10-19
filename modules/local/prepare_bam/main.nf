process PREPARE_BAM {
    label 'process_medium'
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.dict")
    path("tmp/ref_genome.fasta.fai")

    output:
    tuple val(sample_id), path("mapping.bam"),
        path("mapping.bam.csi"),
        path("mapping.raw.flagstat"), emit: bam

    script:
    def rmdup_mode = CLIParamsValidation.rmdup_mode(params.rmdup)
    def rmdup = rmdup_mode == 'none' ? '' : "--rmdup-mode ${rmdup_mode}"
    def left_align = params.left_align_indels.enabled ? '--left-align' : ''
    def platform = params.mapping.platform ?: 'ILLUMINA'
    def library = params.mapping.library ?: sample_id
    def java_options = "-Djava.io.tmpdir=\$TMPDIR " +
        ResourceUtils.javaHeapOption(task.memory)
    def memory_mb = task.memory ? task.memory.toMega() : 0
    def memory_option = memory_mb > 0 ? "--memory-mb ${memory_mb}" : ''
    """
    set -o pipefail

    # Record the input BAM before MAPQ filtering and duplicate handling, so
    # the raw report stays comparable with the mapping path.
    samtools flagstat --threads ${task.cpus} \
        tmp/read1 > mapping.raw.flagstat

    check_bam_reference.py tmp/read1 tmp/ref_genome.fasta.fai \
        --dict tmp/ref_genome.dict

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
        ${memory_option} \
        --java-options "${java_options}" \
        ${rmdup} \
        ${left_align}
    """
}
