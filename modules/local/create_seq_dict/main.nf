process CREATE_SEQ_DICT {
    maxForks 1
    cpus 1
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'

    input:
    path(ref_genome), name: "tmp/ref_genome.fasta"

    output:
    path("ref_genome.dict"), emit: gen_dict

    script:
    """
    gatk CreateSequenceDictionary -R tmp/ref_genome.fasta -O ref_genome.dict
    """
}


