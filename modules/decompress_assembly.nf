process DECOMPRESS_ASSEMBLY {
    maxForks 1
    cpus 1
    afterScript 'stage_cleanup.sh'

    input:
    path(ref_genome), name: "tmp/ref_genome.fasta.gz"

    output:
    path("ref_genome.fasta"), emit: ref_genome

    script:
    """
    bgzip --threads ${task.cpus} --decompress -o ref_genome.fasta tmp/ref_genome.fasta.gz

    """

    stub:
    """
    touch ref_genome.fasta
    """
}