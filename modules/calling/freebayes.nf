process CALL_FREEBAYES {
    cpus params.calling.cpu
    maxForks 1
    errorStrategy 'ignore'
    //afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("freebayes.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    freebayes
    EOF

    parallel -j ${task.cpus} \
    'freebayes --fasta-reference tmp/ref_genome.fasta --targets {} --dont-left-align-indels --ploidy "${params.ploidy}" --min-base-quality "${params.calling.min_base_quality}" ${params.calling.freebayes.extra_args} --bam tmp/input.bam | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools view -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o freebayes.bcf
    """
}