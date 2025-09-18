process CREATE_FAI_INDEX {
    maxForks 1
    cpus 1

    input:
    path(ref_genome)

    output:
    tuple path("reference.fasta"), path("reference.fasta.fai"), emit: fai_index

    script:
    def refName = ref_genome.getName()
    def isGzipped = refName.endsWith('.gz')
    
    """
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Process started - Creating FASTA index"
    
    if [[ "$isGzipped" == "true" ]]; then
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Detected gzipped reference genome, creating ungzipped version"
        gunzip -c ${ref_genome} > reference.fasta
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Ungzipped version created"
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Creating FASTA index"
        samtools faidx --threads ${task.cpus} reference.fasta
    else
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Reference genome is uncompressed"
        mv ${ref_genome} reference.fasta
        echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Creating FASTA index"
        samtools faidx --threads ${task.cpus} reference.fasta
    fi
    
    echo "[\$(date -Iseconds)] [INFO] [CREATE_FAI_INDEX] [${refName}] Process completed - FASTA index created"
    """
}


