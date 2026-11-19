process EXPORT_OUTPUT {
    label 'process_medium'
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    // Only the final export applies output.format. Every file consumed here is
    // an internal BCF produced by GENERATE_CONSENSUS or MERGE_OUTPUTS. The
    // caller passes the target subdirectory so that merged and per-sample
    // exports can run in the same pipeline.
    publishDir(
        { "${params.output.directory}/${publish_subdir}" },
        mode: params.output.publish_mode,
        pattern: '*.bcf',
        enabled: params.output.format == 'bcf'
    )
    publishDir(
        { "${params.output.directory}/${publish_subdir}" },
        mode: params.output.publish_mode,
        pattern: '*.vcf.gz',
        enabled: params.output.format == 'vcf'
    )

    input:
    tuple val(sample_id), val(publish_subdir), path("tmp/input.bcf")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("${sample_id}.bcf"), optional: true, emit: bcf
    tuple val(sample_id), path("${sample_id}.vcf.gz"), optional: true, emit: vcf

    script:
    def extension = params.output.format == 'vcf' ? 'vcf.gz' : 'bcf'
    def output_name = "${sample_id}.${extension}"
    // The per-sample layout publishes under per_sample/ and the merged layout
    // at the output root, so the publish subdirectory also selects which
    // scoped override (output.per_sample or output.merged) applies.
    def scope = publish_subdir == 'per_sample' ? 'per_sample' : 'merged'
    def record_flags = CLIParamsValidation.output_record_flags(
        params.output, scope
    )
    def remove_invariant = record_flags.remove_invariant
        ? '--remove-invariant'
        : ''
    def split_multiallelic = record_flags.split_multiallelic
        ? '--split-multiallelic'
        : ''
    """
    printf '%s\n' tmp/input.bcf > tmp/export_list.txt

    merge_outputs.sh \
        --file-list tmp/export_list.txt \
        --reference tmp/ref_genome.fasta \
        --output ${output_name} \
        --format ${params.output.format} \
        --cpus ${task.cpus} \
        ${remove_invariant} \
        ${split_multiallelic}
    """
}
