process CREATE_SEQ_DICT {
    maxForks 1

    input:
    path(ref_genome)

    output:
    path("${ref_genome.baseName}.dict"), emit: gen_dict

    script:
    def refName = ref_genome.getName()
    
    """
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Process started - Creating sequence dictionary"
    
    gatk CreateSequenceDictionary -R ${ref_genome}

    rm -f ${ref_genome} 2>/dev/null || true
    
    echo "[\$(date -Iseconds)] [INFO] [CREATE_SEQ_DICT] [${refName}] Process completed - Sequence dictionary created"
    """
}


