process CALL_BCFTOOLS {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("bcftools.bcf"), emit: calling_result

    script:
    """
    samtools index --threads ${task.cpus} --csi tmp/input.bam

    cat <<EOF > tmp/sample_id.txt
    bcftools
    EOF

    bcftools mpileup ${params.calling.bcftools.mpileup.extra_args} -Ou --fasta-ref tmp/ref_genome.fasta --threads ${task.cpus} --min-BQ ${params.calling.min_base_quality} --regions-file tmp/coverage.bed tmp/input.bam | bcftools call ${params.calling.bcftools.call.extra_args} -Ou --threads ${task.cpus} | bcftools reheader -s tmp/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter --threads ${task.cpus} -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools view --threads ${task.cpus} -Ou --min-alleles 2 --max-alleles 2 | bcftools norm --threads ${task.cpus} --check-ref w -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o bcftools.bcf
    """
}