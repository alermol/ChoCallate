process CALL_BCFTOOLS {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy 'ignore'
    //afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("bcftools.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/
    
    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    bcftools
    EOF

    parallel -j ${task.cpus} \
    'bcftools mpileup ${params.calling.bcftools.mpileup.extra_args} -Ou --fasta-ref tmp/ref_genome.fasta --threads 1 --min-BQ ${params.calling.min_base_quality} --regions-file {} tmp/input.bam | bcftools call ${params.calling.bcftools.call.extra_args} -Ou --threads 1 | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt | bcftools norm --check-ref w -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools view -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o bcftools.bcf
    """
}