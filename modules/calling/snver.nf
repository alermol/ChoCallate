process CALL_SNVER {
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

    output:
    tuple val(sample_id), path("snver.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    snver
    EOF

    parallel -j ${task.cpus} \
    "mkdir -p tmp/calling_chunks/{#}; snver ${params.calling.snver.extra_args} -i tmp/input.bam -r tmp/ref_genome.fasta -l {} -bq "${params.calling.min_base_quality}" -n "${params.ploidy}" -o tmp/calling_chunks/{#}/{#}
    bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai tmp/calling_chunks/{#}/{#}.filter.vcf | bgzip --threads 1 -c > tmp/calling_chunks/{#}/{#}.filter.vcf.gz
    tabix --threads 1 --csi tmp/calling_chunks/{#}/{#}.filter.vcf.gz
    bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai tmp/calling_chunks/{#}/{#}.indel.filter.vcf | bgzip --threads 1 -c > tmp/calling_chunks/{#}/{#}.indel.filter.vcf.gz
    tabix --threads 1 --csi tmp/calling_chunks/{#}/{#}.indel.filter.vcf.gz
    bcftools concat --allow-overlaps --threads 1 -Ou tmp/calling_chunks/{#}/{#}.filter.vcf.gz tmp/calling_chunks/{#}/{#}.indel.filter.vcf.gz | bcftools sort -Ou | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}/{#}.bcf" ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*/*.bcf | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools view -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o snver.bcf
    """
}