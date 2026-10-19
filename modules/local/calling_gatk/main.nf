process CALLING_GATK {
    label 'process_high'
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam"), path("tmp/input.bam.csi"), path("tmp/bed_chunks/*", arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/ref_genome.dict")

    output:
    tuple val(sample_id), path("gatk.bcf"), emit: calling_result

    script:
    def caller = params.calling.gatk ?: [:]
    def extra_args = caller.extra_args ?: ''
    def interval_padding = caller.interval_padding != null ? caller.interval_padding : 150
    def min_snp = caller.min_snp_qual
    def min_indel = caller.min_indel_qual
    def heap_option = ResourceUtils.javaHeapOption(task.memory, task.cpus)
    assert min_snp != null : 'calling.gatk.min_snp_qual is required'
    assert min_indel != null : 'calling.gatk.min_indel_qual is required'
    def call_confidence = [min_snp, min_indel].min()
    def qual_filter = "(TYPE=\\\"snp\\\" && QUAL<${min_snp}) || (TYPE!=\\\"snp\\\" && QUAL<${min_indel})"
    """
    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    # The padding lets GATK assemble across chunk edges, but it also makes
    # neighbouring chunks call the same boundary sites (possibly with
    # different allele representations, which -d exact cannot deduplicate).
    # Restrict each chunk to its own unpadded interval before concatenating;
    # -T streams the BED and needs no index.
    cat > tmp/run_chunk.sh <<'EOF'
    set -euo pipefail
    gatk HaplotypeCaller \
        --standard-min-confidence-threshold-for-calling ${call_confidence} \
        ${extra_args} \
        --java-options "-Djava.io.tmpdir=\$TMPDIR ${heap_option}" \
        --do-not-run-physical-phasing true \
        --smith-waterman FASTEST_AVAILABLE \
        --create-output-variant-index false \
        --native-pair-hmm-threads 1 \
        --input tmp/input.bam --reference tmp/ref_genome.fasta \
        --output /dev/stdout --intervals "\$1" \
        --interval-padding ${interval_padding} \
        --sequence-dictionary tmp/ref_genome.dict \
        --tmp-dir "\$TMPDIR" \
    | bcftools reheader -s tmp/calling_chunks/sample_id.txt \
        -f tmp/ref_genome.fasta.fai \
    | bcftools view -T "\$1" -Ob -o "tmp/calling_chunks/\$2.bcf"
    EOF

    # SIGKILLed children appear as signal 9, exit 137, or GATK's 247.
    parallel --halt now,fail=1 --joblog tmp/parallel.log \
        -j ${task.cpus} bash tmp/run_chunk.sh {} {#} \
        ::: tmp/bed_chunks/*.bed || {
        status=\$?
        if awk 'NR > 1 && (\$7 == 137 || \$7 == 247 || \$8 == 9) {
            found = 1 } END { exit !found }' tmp/parallel.log; then
            exit 137
        fi
        exit "\${status}"
    }

    # Each chunk now contributes only its own interval, so the concatenated
    # records are disjoint and ordered. Sort only this caller as a guard for
    # any remaining boundary order; the exact chunk targets of
    # bcftools/freebayes are already ordered.
    bcftools concat --naive --threads ${task.cpus} -Ob \$(ls -v tmp/calling_chunks/*.bcf) \
    | bcftools norm --threads ${task.cpus} -f tmp/ref_genome.fasta -m -any --atomize --check-ref x -Ou \
    | bcftools sort -Ou \
    | bcftools norm --threads ${task.cpus} -d exact -Ou \
    | bcftools filter -s LowQual -m + -e"${qual_filter}" -Ob -o gatk.bcf
    """
}
