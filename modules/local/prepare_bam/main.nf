process PREPARE_BAM {
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
    tuple val(sample_id), path("output.bam"), emit: bam
    
    script:
    def rmdup = params.rmdup.enabled ? '| picard MarkDuplicates -I - -O - -M /dev/null --VALIDATION_STRINGENCY SILENT --REMOVE_DUPLICATES true' : ''
    def lai = params.left_align_indels.enabled ? '| gatk LeftAlignIndels -I /dev/stdin -O /dev/stdout -R tmp/ref_genome.fasta --sequence-dictionary tmp/ref_genome.dict -OBI false' : ''
    def fixmate = params.input.reads_type != 'se' ? "| samtools sort -n - | samtools fixmate --threads ${task.cpus} -m - -" : ''
    def reads = "tmp/read1"
    """
    samtools view -F 4 --threads ${task.cpus} -b -q ${params.bam_filter.min_map_qual} ${reads} ${fixmate} | \
        samtools sort --threads ${task.cpus} ${rmdup} ${lai} | \
        samtools addreplacerg -w -r "@RG\tID:${sample_id}\tSM:${sample_id}" --threads ${task.cpus} -o mapping.bam
    """
}