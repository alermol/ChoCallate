process MAPPING_MINIMAP2 {
    maxForks 1
    cpus params.mapping.cpu
    beforeScript 'export TMPDIR=$(mktemp -d -p $PWD/)'
    afterScript 'stage_cleanup.sh'
    
    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/read1"), path("tmp/read2"), path("tmp/read3")
    path("tmp/ref_genome_real.fasta")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.dict")
    path("tmp/ref_genome.fasta.fai")

    output:
    tuple val(sample_id), path("mapping.bam"), emit: bam

    script:
    def rg_id = "@RG\\tID:${sample_id}\\tSM:${sample_id}"
    def rmdup = params.rmdup.enabled ? 'picard MarkDuplicates -I tmp/stage1.bam -O tmp/stage2.bam -M /dev/null --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true' : 'mv tmp/stage1.bam tmp/stage2.bam'
    def lai = params.left_align_indels.enabled ? 'gatk LeftAlignIndels -I tmp/stage2.bam -O mapping.bam -R tmp/ref_genome.fasta --sequence-dictionary tmp/ref_genome.dict -OBI false' : 'mv tmp/stage2.bam mapping.bam'
    def fixmate = params.input.reads_type != 'se' ? "| samtools sort --threads ${task.cpus} -n | samtools fixmate --threads ${task.cpus} -m - -" : ''
    def reads = params.input.reads_type == 'se' ? "tmp/read1" : 
                params.input.reads_type == 'pe' ? "tmp/read1 tmp/read2" : ''
    """
    GENOME_BASENAME=\$(basename \$(realpath tmp/ref_genome.fasta))
    
    minimap2 \
        ${params.mapping.extra_args} 
        -ax sr \
        -t ${task.cpus} \
        -R "${rg_id}" 
        -a ${params.input.reference_index_dir}/\${GENOME_BASENAME}.mmi \
        ${reads} | \
    samtools view -F 4 --threads ${task.cpus} -b -q ${params.bam_filter.min_map_qual} ${fixmate} | \
    samtools sort --threads ${task.cpus} -o tmp/stage1.bam
    
    ${rmdup}
     
    ${lai}
    """
}