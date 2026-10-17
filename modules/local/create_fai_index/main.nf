process CREATE_FAI_INDEX {
    cpus Runtime.runtime.availableProcessors()
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    input:
    path(ref_genome), name: "tmp/ref_genome.fasta"

    output:
    path("ref_genome.fasta.fai"), emit: fai_index

    script:
    """
    samtools faidx --threads ${task.cpus} --output ref_genome.fasta.fai tmp/ref_genome.fasta
    """
}

