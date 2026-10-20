process MERGE_OUTPUTS {
    label 'process_medium'
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "consensus"

    input:
    // The sample ids are carried explicitly so the merge order follows
    // samples_tsv even though `dir/*` staging preserves the file names.
    tuple val(sample_ids), path('tmp/bcfs/*', arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    val dropped_samples

    output:
    // Internal merged representation.  Filtering, re-joining and the final
    // BCF/VCF conversion are done by EXPORT_OUTPUT.
    path("consensus.bcf"), emit: consensus

    script:
    def sample_id_list = sample_ids.split(',').collect { id -> id.trim() }
    def merge_paths = sample_id_list.collect { id -> "tmp/bcfs/${id}.bcf" }.join(' ')
    def dropped = dropped_samples ? dropped_samples.join(',') : ''
    def dropped_arg = dropped
        ? "--samples-dropped '${dropped}'"
        : ''
    """
    printf '%s\n' ${merge_paths} > tmp/merge_list.txt

    # Keep one record per ALT in the internal merge.  EXPORT_OUTPUT re-joins
    # multi-allelic records when output.split_multiallelic is false.
    merge_outputs.sh \
        --file-list tmp/merge_list.txt \
        --reference tmp/ref_genome.fasta \
        --output consensus.bcf \
        --format bcf \
        --split-multiallelic \
        --cpus ${task.cpus} \
        ${dropped_arg}
    """
}
