process GENERATE_CONSENSUS {
    label 'process_medium'
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    // The `dir/*` staging pattern preserves the caller file names (bcftools.bcf,
    // freebayes.bcf, gatk.bcf), so the consensus step knows which caller each
    // file came from even when a caller is missing.
    tuple val(sample_id), val(received_callers), path('tmp/callers/*', arity: '1..*'), path("tmp/input.bam"), path("tmp/input.bam.csi"), path("tmp/bed_chunks/*", arity: '1..*')
    val configured_callers
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")

    output:
    // Internal representation only. The requested output format is applied by
    // EXPORT_OUTPUT after the optional multi-sample merge.
    tuple val(sample_id), path("${sample_id}.bcf"), emit: consensus_bcf


    script:
    def configured = configured_callers.split(',').collect { caller -> caller.trim() }.findAll { caller -> caller }
    def received = received_callers.split(',').collect { caller -> caller.trim() }.findAll { caller -> caller }
    def failed = configured.findAll { caller -> !(caller in received) }
    if (failed && params.consensus.threshold > received.size()) {
        error "Sample ${sample_id}: consensus.threshold (${params.consensus.threshold}) exceeds the ${received.size()} surviving caller(s) [${received.join(',')}] after failure(s) of ${failed.join(',')}"
    }
    def input_args = received.collect { caller -> "--input ${caller}=tmp/callers/${caller}.bcf" }.join(' ')
    def compress_blocks = params.output.compress_reference_blocks
        ? "--compress_reference_blocks"
        : ""
    // Per-sample export discards reference records when remove_invariant is
    // true, so skip writing them at all.  A merged output still needs them for
    // bcftools merge --gvcf to assign ALT depths across samples.
    def variants_only = params.output.type == 'sample'
        && params.output.remove_invariant
    def variants_only_arg = variants_only ? "--variants_only" : ""
    def absent_vote = params.consensus.absent_vote != null ? params.consensus.absent_vote : 'ref'
    def filtered_vote = params.consensus.filtered_vote != null
        ? params.consensus.filtered_vote
        : 'abstain'
    def tie_policy = params.consensus.tie_policy != null ? params.consensus.tie_policy : 'nocall'
    def min_coverage = params.coverage.min_coverage != null
        ? params.coverage.min_coverage
        : 5
    """
    mkdir -p tmp/bcf_chunks/

    parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' \
        ::: tmp/callers/*.bcf

    cat > tmp/run_chunk.sh <<'EOF'
    bgzip --threads 1 "\$1"
    tabix --threads 1 --csi -p bed "\$1.gz"
    generate_consensus.py ${input_args} \
        --callers "${configured_callers}" \
        --bed "\$1.gz" \
        --bam tmp/input.bam \
        --output "tmp/bcf_chunks/\$2.bcf" \
        --sample_name "${sample_id}" \
        --reference tmp/ref_genome.fasta \
        --consensus_threshold "${params.consensus.threshold}" \
        --absent_vote "${absent_vote}" \
        --filtered_vote "${filtered_vote}" \
        --tie_policy "${tie_policy}" \
        --min_coverage "${min_coverage}" \
        --version "${workflow.manifest.version}" \
        ${compress_blocks} \
        ${variants_only_arg}
    EOF

    parallel -j ${task.cpus} bash tmp/run_chunk.sh {} {#} \
        ::: tmp/bed_chunks/*.bed

    # The interval chunks are consumed in BED order; ls -v preserves that order
    # through concat, so no full sort is needed afterwards.  The per-sample
    # consensus stays split by ALT; reference records are written unless the
    # per-sample export would drop them.  Output filtering and re-joining happen
    # once in EXPORT_OUTPUT.
    bcftools concat --naive --threads ${task.cpus} -Ob \
        -o ${sample_id}.bcf \$(ls -v tmp/bcf_chunks/*.bcf)
    """
}
