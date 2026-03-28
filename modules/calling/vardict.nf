process CALL_VARDICT {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy 'ignore'
    //afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("vardict.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    vardict
    EOF

    parallel -j ${task.cpus} \
    'vardict-java ${params.calling.vardict.extra_args} -G tmp/ref_genome.fasta -b tmp/input.bam -fisher -th 1 -q "${params.calling.min_base_quality}" -VS SILENT -c 1 -S 2 -E 3 -g 4 {} | var2vcf_valid.pl -S -q "${params.calling.min_base_quality}" -E | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools view -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o vardict.bcf
    """
}