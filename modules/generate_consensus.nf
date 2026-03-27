process GENERATE_CONSENSUS {
    maxForks 1
    cpus params.consensus.cpu
    //afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    publishDir "${params.output.directory}/per_sample", mode: 'move', pattern: '*.bcf', enabled: params.output.type == 'sample' && params.output.format == 'bcf'

    input:
    tuple val(sample_id), path('tmp/?.bcf', arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("${sample_id}.bcf"), emit: consensus

    script:
    """
    mkdir -p tmp/bed_chunks/
    mkdir -p tmp/vcf_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf

    parallel -j ${task.cpus} \
    'bgzip --threads 1 {}
    tabix --threads 1 --csi -p bed {}.gz
    generate_consensus.py --input tmp/*.bcf --bed {}.gz --output tmp/vcf_chunks/{#}.bcf --sample_name "${sample_id}" --reference tmp/ref_genome.fasta --consensus_threshold "${params.consensus.threshold}"
    bcftools index --threads 1 --csi tmp/vcf_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/vcf_chunks/*.bcf | bcftools sort -Ob -o "${sample_id}.bcf"
    """
}


