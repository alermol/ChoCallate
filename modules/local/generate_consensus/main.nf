process GENERATE_CONSENSUS {
    cpus params.consensus.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    publishDir "${params.output.directory}/per_sample", mode: 'move', pattern: '*.bcf', enabled: params.output.type == 'sample' && params.output.format == 'bcf'
    publishDir "${params.output.directory}/per_sample", mode: 'move', pattern: '*.vcf.gz', enabled: params.output.type == 'sample' && params.output.format == 'vcf'

    input:
    tuple val(sample_id), path('tmp/?.bcf', arity: '1..*')
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("${sample_id}.bcf"), emit: consensus_bcf, optional: true
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: consensus_vcf, optional: true
    

    script:
    def output_format = params.output.format == 'vcf' ? "-Oz -o ${sample_id}.vcf.gz" : "-Ob -o ${sample_id}.bcf"
    def filter_invariant = params.output.remove_invariant && params.output.type == 'sample' ? "| bcftools filter --threads ${task.cpus} -e 'COUNT(GT=\"hom\")=N_SAMPLES || COUNT(GT=\"het\")=N_SAMPLES' -Ou" : ""
    def split_multiallelic = params.output.split_multiallelic ? "--split_multiallelic" : ""
    def remove_invariant = params.output.remove_invariant ? "--remove_invariant" : ""
    """
    mkdir -p tmp/bed_chunks/
    mkdir -p tmp/vcf_chunks/
    bedops --chop 10000 tmp/coverage.bed > tmp/coverage.bed.chopped
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed.chopped tmp/bed_chunks/

    parallel -j ${task.cpus} 'bcftools index --threads 1 --csi {}' ::: tmp/*.bcf

    parallel -j ${task.cpus} \
    'bgzip --threads 1 {}
    tabix --threads 1 --csi -p bed {}.gz
    generate_consensus.py --input tmp/*.bcf --bed {}.gz --output tmp/vcf_chunks/{#}.bcf --sample_name "${sample_id}" --reference tmp/ref_genome.fasta --consensus_threshold "${params.consensus.threshold}" --version "${workflow.manifest.version}" ${split_multiallelic} ${remove_invariant}
    bcftools index --threads 1 --csi tmp/vcf_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --naive --threads ${task.cpus} -Ob tmp/vcf_chunks/*.bcf ${filter_invariant} | bcftools sort ${output_format}
    """
}


