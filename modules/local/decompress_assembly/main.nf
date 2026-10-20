process DECOMPRESS_ASSEMBLY {
    label 'process_low'
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    input:
    path(ref_genome), name: "tmp/ref_genome.fasta.gz"

    output:
    path("ref_genome.fasta"), emit: ref_genome

    script:
    """
    bgzip --threads ${task.cpus} --decompress \
        -o ref_genome.fasta tmp/ref_genome.fasta.gz
    """
}