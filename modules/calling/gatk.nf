process CALL_GATK {
    cpus params.calling.cpu
    maxForks 1
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")
    path("tmp/ref_genome.dict")

    output:
    tuple val(sample_id), path("gatk.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'gatk HaplotypeCaller ${params.calling.gatk.extra_args} --do-not-run-physical-phasing true --smith-waterman FASTEST_AVAILABLE --create-output-variant-index false --native-pair-hmm-threads 1 --input tmp/input.bam --reference tmp/ref_genome.fasta --output /dev/stdout --intervals {} --sequence-dictionary tmp/ref_genome.dict | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools norm --check-ref x -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools view --threads ${task.cpus} -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o gatk.bcf
    """
}