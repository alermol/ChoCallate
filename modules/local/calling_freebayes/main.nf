process CALLING_FREEBAYES {
    label 'process_medium'
    cpus params.calling.cpu
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam"), path("tmp/input.bam.csi"), path("tmp/bed_chunks/*", arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/ref_genome.dict")

    output:
    tuple val(sample_id), path("freebayes.bcf"), emit: calling_result

    script:
    def caller = params.calling.freebayes ?: [:]
    def extra_args = caller.extra_args ?: ''
    def use_best_n_alleles = caller.use_best_n_alleles
    def limit_coverage = caller.limit_coverage
    def min_snp = caller.min_snp_qual
    def min_indel = caller.min_indel_qual
    assert use_best_n_alleles != null :
        'calling.freebayes.use_best_n_alleles is required'
    assert limit_coverage != null :
        'calling.freebayes.limit_coverage is required'
    assert min_snp != null : 'calling.freebayes.min_snp_qual is required'
    assert min_indel != null : 'calling.freebayes.min_indel_qual is required'
    def qual_filter = "(TYPE=\\\"snp\\\" && QUAL<${min_snp}) || (TYPE!=\\\"snp\\\" && QUAL<${min_indel})"
    """
    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    cat > tmp/run_chunk.sh <<'EOF'
    set -euo pipefail
    freebayes --fasta-reference tmp/ref_genome.fasta --targets "\$1" \
        --use-best-n-alleles ${use_best_n_alleles} \
        --limit-coverage ${limit_coverage} \
        --hwe-priors-off --no-population-priors \
        --dont-left-align-indels ${extra_args} --bam tmp/input.bam \
    | bcftools reheader -s tmp/calling_chunks/sample_id.txt \
        -f tmp/ref_genome.fasta.fai \
    | bcftools view -Ob -o "tmp/calling_chunks/\$2.bcf"
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

    # The chunks are emitted in genomic order; ls -v keeps that order for
    # concat, so no full sort is needed afterwards.
    bcftools concat --naive --threads ${task.cpus} -Ob \$(ls -v tmp/calling_chunks/*.bcf) \
    | bcftools norm --threads ${task.cpus} -f tmp/ref_genome.fasta -m -any --atomize --check-ref x -Ou \
    | bcftools norm --threads ${task.cpus} -d exact -Ou \
    | bcftools filter -s LowQual -m + -e"${qual_filter}" -Ob -o freebayes.bcf
    """
}
