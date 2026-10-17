process MERGE_OUTPUTS {
    label 'process_medium'
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    publishDir(
        "${params.output.directory}",
        mode: params.output.publish_mode,
        pattern: 'consensus.bcf',
        enabled: params.output.type == 'single' &&
                 params.output.format == 'bcf'
    )
    publishDir(
        "${params.output.directory}",
        mode: params.output.publish_mode,
        pattern: 'consensus.vcf.gz',
        enabled: params.output.type == 'single' &&
                 params.output.format == 'vcf'
    )

    input:
    // The sample ids are carried explicitly so the merge order follows
    // samples_tsv even though `dir/*` staging preserves the file names.
    tuple val(sample_ids), path('tmp/bcfs/*', arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")

    output:
    path("consensus.bcf"), optional: true
    path("consensus.vcf.gz"), optional: true

    script:
    def sample_id_list = sample_ids.split(',').collect { id -> id.trim() }
    def merge_paths = sample_id_list.collect { id -> "tmp/bcfs/${id}.bcf" }.join(' ')
    def output_name = params.output.format == 'vcf' ? 'consensus.vcf.gz' : 'consensus.bcf'
    def remove_invariant = params.output.remove_invariant ? '--remove-invariant' : ''
    def split_multiallelic = params.output.split_multiallelic ? '--split-multiallelic' : ''
    """
    printf '%s\n' ${merge_paths} > tmp/merge_list.txt

    merge_outputs.sh \
        --file-list tmp/merge_list.txt \
        --reference tmp/ref_genome.fasta \
        --output ${output_name} \
        --format ${params.output.format} \
        --cpus ${task.cpus} \
        ${remove_invariant} \
        ${split_multiallelic}
    """
}
