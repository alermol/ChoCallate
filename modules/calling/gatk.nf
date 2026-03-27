process CALL_GATK {
    cpus params.calling.cpu
    maxForks 1
    errorStrategy 'ignore'
    //afterScript 'stage_cleanup.sh'

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
    'gatk HaplotypeCaller ${params.calling.gatk.extra_args} --input tmp/input.bam --reference tmp/ref_genome.fasta --min-base-quality-score "${params.calling.min_base_quality}" --output /dev/stdout --intervals {} --ploidy "${params.ploidy}" --sequence-dictionary tmp/ref_genome.dict | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools view -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o gatk.bcf
    """
}