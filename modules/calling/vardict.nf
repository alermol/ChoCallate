process CALL_VARDICT {
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
    tuple val(sample_id), path("vardict.bcf"), emit: calling_result

    script:
    """
    samtools index --threads ${task.cpus} --csi tmp/input.bam

    cat <<EOF > tmp/sample_id.txt
    vardict
    EOF

    vardict-java ${params.calling.vardict.extra_args} --nosv -k 0 -G tmp/ref_genome.fasta -b tmp/input.bam -fisher -th ${task.cpus} -VS SILENT -c 1 -S 2 -E 3 -g 4 -th ${task.cpus} tmp/coverage.bed | var2vcf_valid.pl -S -E | bcftools reheader -s tmp/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools view -Ou --min-alleles 2 --max-alleles 2 | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o vardict.bcf
    """
}