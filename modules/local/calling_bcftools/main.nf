process CALLING_BCFTOOLS {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy 'ignore'
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
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
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'bcftools mpileup ${params.calling.bcftools.mpileup.extra_args} --count-orphans -Ou --fasta-ref tmp/ref_genome.fasta --regions-file {} tmp/input.bam | bcftools call ${params.calling.bcftools.call.extra_args} --multiallelic-caller --variants-only -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools norm --check-ref x -Ob --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --naive --threads ${task.cpus} -Ob tmp/calling_chunks/*.bcf | bcftools view --threads ${task.cpus} -Ou --min-alleles 2 --max-alleles 2 | bcftools sort -Ou -o bcftools.bcf
    """
}