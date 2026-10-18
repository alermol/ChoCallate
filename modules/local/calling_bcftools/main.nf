process CALLING_BCFTOOLS {
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
    tuple val(sample_id), path("bcftools.bcf"), emit: calling_result

    script:
    def caller = params.calling.bcftools ?: [:]
    def call_args = caller.call?.extra_args ?: ''
    def mpileup_args = caller.mpileup?.extra_args ?: ''
    def min_snp = caller.min_snp_qual
    def min_indel = caller.min_indel_qual
    assert min_snp != null : 'calling.bcftools.min_snp_qual is required'
    assert min_indel != null : 'calling.bcftools.min_indel_qual is required'
    def qual_filter = "(TYPE=\\\"snp\\\" && QUAL<${min_snp}) || (TYPE!=\\\"snp\\\" && QUAL<${min_indel})"
    """
    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'bcftools mpileup ${mpileup_args} -a FORMAT/AD,FORMAT/DP --count-orphans -Ou --fasta-ref tmp/ref_genome.fasta --regions-file {} tmp/input.bam | bcftools call ${call_args} --multiallelic-caller --variants-only -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -s LowQual -m + -e"${qual_filter}" -Ob -o tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    # The chunks are emitted in genomic order; ls -v keeps that order for
    # concat, so no full sort is needed afterwards.
    bcftools concat --naive --threads ${task.cpus} -Ob \$(ls -v tmp/calling_chunks/*.bcf) \
    | bcftools norm --threads ${task.cpus} -f tmp/ref_genome.fasta -m -any --atomize --check-ref x -Ou \
    | bcftools norm --threads ${task.cpus} -d exact -Ob -o bcftools.bcf
    """
}
