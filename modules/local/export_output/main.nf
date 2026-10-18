process EXPORT_OUTPUT {
    label 'process_medium'
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    // Only the final export uses output.format and output.type. Every file
    // consumed here is an internal BCF produced by GENERATE_CONSENSUS or
    // MERGE_OUTPUTS.
    publishDir(
        params.output.type == 'single'
            ? params.output.directory
            : "${params.output.directory}/per_sample",
        mode: params.output.publish_mode,
        pattern: '*.bcf',
        enabled: params.output.format == 'bcf'
    )
    publishDir(
        params.output.type == 'single'
            ? params.output.directory
            : "${params.output.directory}/per_sample",
        mode: params.output.publish_mode,
        pattern: '*.vcf.gz',
        enabled: params.output.format == 'vcf'
    )

    input:
    tuple val(sample_id), path("tmp/input.bcf")
    path("tmp/ref_genome.fasta")

    output:
    tuple val(sample_id), path("${sample_id}.bcf"), optional: true, emit: bcf
    tuple val(sample_id), path("${sample_id}.vcf.gz"), optional: true, emit: vcf

    script:
    def extension = params.output.format == 'vcf' ? 'vcf.gz' : 'bcf'
    def output_name = "${sample_id}.${extension}"
    def remove_invariant = params.output.remove_invariant
        ? '--remove-invariant'
        : ''
    def split_multiallelic = params.output.split_multiallelic
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
