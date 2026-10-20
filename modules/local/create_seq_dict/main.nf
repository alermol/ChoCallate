process CREATE_SEQ_DICT {
    label 'process_low'
    cpus 1
    beforeScript 'export TMPDIR="$(mktemp -d -p "${PWD}")"'
    afterScript 'stage_cleanup.sh "${nxf_main_ret:-0}"'

    input:
    path(ref_genome), name: "tmp/ref_genome.fasta"

    output:
    path("ref_genome.dict"), emit: gen_dict

    script:
    def heap_option = ResourceUtils.javaHeapOption(task.memory)
    """
    gatk CreateSequenceDictionary \
        --java-options "-Djava.io.tmpdir=\$TMPDIR ${heap_option}" \
        --TMP_DIR "\$TMPDIR" \
        -R tmp/ref_genome.fasta -O ref_genome.dict
    """
}


