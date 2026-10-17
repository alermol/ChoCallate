process CALLING_GATK {
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
    tuple val(sample_id), path("gatk.bcf"), emit: calling_result

    script:
    def caller = params.calling.gatk ?: [:]
    def extra_args = caller.extra_args ?: ''
    def interval_padding = caller.interval_padding != null ? caller.interval_padding : 150
    def min_snp = caller.min_snp_qual
    def min_indel = caller.min_indel_qual
    assert min_snp != null : 'calling.gatk.min_snp_qual is required'
    assert min_indel != null : 'calling.gatk.min_indel_qual is required'
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
    'gatk HaplotypeCaller ${extra_args} --java-options "-Djava.io.tmpdir=\$TMPDIR" --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE --create-output-variant-index false --native-pair-hmm-threads 1 --input tmp/input.bam --reference tmp/ref_genome.fasta --output /dev/stdout --intervals {} --interval-padding ${interval_padding} --sequence-dictionary tmp/ref_genome.dict --tmp-dir \$TMPDIR | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"${qual_filter}" -Ob -o tmp/calling_chunks/{#}.bcf
    bcftools index --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    # Split multi-allelics and atomise complex alleles once, before voting, so
    # every caller uses the same representation.  The sort + exact dedup remove
    # records duplicated at chunk boundaries.
    bcftools concat --naive --threads ${task.cpus} -Ob tmp/calling_chunks/*.bcf \
    | bcftools norm --threads ${task.cpus} -f tmp/ref_genome.fasta -m -any --atomize --check-ref x -Ou \
    | bcftools sort -Ou \
    | bcftools norm --threads ${task.cpus} -d exact -Ob -o gatk.bcf
    """
}
