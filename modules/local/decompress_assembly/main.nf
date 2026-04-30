process DECOMPRESS_ASSEMBLY {
    cpus Runtime.runtime.availableProcessors()
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    input:
    path(ref_genome), name: "tmp/ref_genome.fasta.gz"

    output:
    path("ref_genome.fasta"), emit: ref_genome

    script:
    """
    bgzip --threads ${task.cpus} --decompress -o ref_genome.fasta tmp/ref_genome.fasta.gz
    """
}