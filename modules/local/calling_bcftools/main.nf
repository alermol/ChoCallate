process CALLING_BCFTOOLS {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy { params.calling.allow_caller_failure ? 'ignore' : (task.attempt <= 2 ? 'retry' : 'terminate') }
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam"), path("tmp/coverage.bed")
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
    mkdir -p tmp/bed_chunks/
    mkdir -p tmp/calling_chunks/

    gatk SplitIntervals --java-options "-Djava.io.tmpdir=\$TMPDIR" --tmp-dir \$TMPDIR -R tmp/ref_genome.fasta -L tmp/coverage.bed --scatter-count ${task.cpus} -O tmp/bed_chunks/
    parallel -j ${task.cpus} 'gatk IntervalListToBed -I {} -O {//}/{/.}.bed.tmp; cut -f 1-3 {//}/{/.}.bed.tmp | tee {//}/{/.}.bed; rm {} {//}/{/.}.bed.tmp' ::: tmp/bed_chunks/*

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'bcftools mpileup ${mpileup_args} -a FORMAT/AD,FORMAT/DP --count-orphans -Ou --fasta-ref tmp/ref_genome.fasta --regions-file {} tmp/input.bam | bcftools call ${call_args} --multiallelic-caller --variants-only -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"${qual_filter}" -Ou | bcftools norm --check-ref x -Ob --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf;
    bcftools index --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --naive --threads ${task.cpus} -Ob tmp/calling_chunks/*.bcf | bcftools view --threads ${task.cpus} -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o bcftools.bcf
    """
}
