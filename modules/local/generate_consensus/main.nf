process GENERATE_CONSENSUS {
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    publishDir "${params.output.directory}/per_sample", mode: 'move', pattern: '*.bcf', enabled: params.output.type == 'sample' && params.output.format == 'bcf'
    publishDir "${params.output.directory}/per_sample", mode: 'move', pattern: '*.vcf.gz', enabled: params.output.type == 'sample' && params.output.format == 'vcf'

    input:
    tuple val(sample_id), path('tmp/?.bcf', arity: '1..*'), path("tmp/input.bam"), path("tmp/coverage.bed")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/ref_genome.dict")

    output:
    tuple val(sample_id), path("${sample_id}.bcf"), emit: consensus_bcf, optional: true
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: consensus_vcf, optional: true


    script:
    def output_format = params.output.format == 'vcf' ? "-Oz -o ${sample_id}.vcf.gz" : "-Ob -o ${sample_id}.bcf"
    def filter_invariant = params.output.remove_invariant && params.output.type == 'sample' ? "| bcftools filter --threads ${task.cpus} -e 'COUNT(GT=\"RR\")=N_SAMPLES' -Ou" : ""
    def split_multiallelic = params.output.split_multiallelic ? "--split_multiallelic" : ""
    def remove_invariant = params.output.remove_invariant ? "--remove_invariant" : ""
    def fill_tags = "AN,AC,AF,NS,AC_Hom,AC_Het,MAF,TYPE,F_MISSING,'DP:1=int(sum(FORMAT/DP))'"
    """
    mkdir -p tmp/bed_chunks/
    mkdir -p tmp/vcf_chunks/

    gatk SplitIntervals --java-options "-Djava.io.tmpdir=\$TMPDIR" --tmp-dir \$TMPDIR -R tmp/ref_genome.fasta -L tmp/coverage.bed --scatter-count ${task.cpus} -O tmp/bed_chunks/
    parallel -j ${task.cpus} 'gatk IntervalListToBed -I {} -O {//}/{/.}.bed.tmp; cut -f 1-3 {//}/{/.}.bed.tmp | tee {//}/{/.}.bed; rm {} {//}/{/.}.bed.tmp' ::: tmp/bed_chunks/*

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf

    parallel -j ${task.cpus} \
    'bgzip --threads 1 {}
    tabix --threads 1 --csi -p bed {}.gz
    generate_consensus.py --input tmp/*.bcf --bed {}.gz --bam tmp/input.bam --output tmp/vcf_chunks/{#}.bcf --sample_name "${sample_id}" --reference tmp/ref_genome.fasta --consensus_threshold "${params.consensus.threshold}" --version "${workflow.manifest.version}" ${split_multiallelic} ${remove_invariant}
    bcftools index --threads 1 --csi tmp/vcf_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --naive --threads ${task.cpus} -Ob tmp/vcf_chunks/*.bcf ${filter_invariant} \
    | bcftools +fill-tags -Ou -- -t ${fill_tags} \
    | bcftools sort ${output_format}
    """
}
