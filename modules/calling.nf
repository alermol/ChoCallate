process CALL_BCFTOOLS {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("bcftools.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/
    
    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'bcftools mpileup ${params.calling.bcftools.mpileup.extra_args} -Ou --fasta-ref tmp/ref_genome.fasta --threads 1 --min-BQ ${params.calling.min_base_quality} --regions-file {} tmp/input.bam | bcftools call ${params.calling.bcftools.call.extra_args} -Ou --threads 1 | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt | bcftools norm --check-ref w -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools sort -Ou -o bcftools.bcf
    """
}

process CALL_FREEBAYES {
    cpus params.calling.cpu
    maxForks 1
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("freebayes.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'freebayes --fasta-reference tmp/ref_genome.fasta --targets {} --dont-left-align-indels --ploidy "${params.ploidy}" --min-base-quality "${params.calling.min_base_quality}" ${params.calling.freebayes.extra_args} --bam tmp/input.bam | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools sort -Ou -o freebayes.bcf
    """
}

process CALL_GATK {
    cpus params.calling.cpu
    maxForks 1
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")
    path("tmp/ref_genome.dict")

    output:
    tuple val(sample_id), path("gatk.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'gatk HaplotypeCaller ${params.calling.gatk.extra_args} --input tmp/input.bam --reference tmp/ref_genome.fasta --min-base-quality-score "${params.calling.min_base_quality}" --output /dev/stdout --intervals {} --ploidy "${params.ploidy}" --sequence-dictionary tmp/ref_genome.dict | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools reheader -s tmp/calling_chunks/sample_id.txt | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools sort -Ou -o gatk.bcf
    """
}

process CALL_SNVER {
    cpus params.calling.cpu
    maxForks 1
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("snver.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    "mkdir -p tmp/calling_chunks/{#}; snver ${params.calling.snver.extra_args} -i tmp/input.bam -r tmp/ref_genome.fasta -l {} -bq "${params.calling.min_base_quality}" -n "${params.ploidy}" -o tmp/calling_chunks/{#}/{#}
    bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai tmp/calling_chunks/{#}/{#}.filter.vcf | bgzip --threads 1 -c > tmp/calling_chunks/{#}/{#}.filter.vcf.gz
    tabix --threads 1 --csi tmp/calling_chunks/{#}/{#}.filter.vcf.gz
    bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai tmp/calling_chunks/{#}/{#}.indel.filter.vcf | bgzip --threads 1 -c > tmp/calling_chunks/{#}/{#}.indel.filter.vcf.gz
    tabix --threads 1 --csi tmp/calling_chunks/{#}/{#}.indel.filter.vcf.gz
    bcftools concat --allow-overlaps --threads 1 -Ou tmp/calling_chunks/{#}/{#}.filter.vcf.gz tmp/calling_chunks/{#}/{#}.indel.filter.vcf.gz | bcftools sort -Ou | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}/{#}.bcf" ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*/*.bcf | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools sort -Ou -o snver.bcf
    """
}

process CALL_VARDICT {
    maxForks 1
    cpus params.calling.cpu
    errorStrategy 'ignore'
    afterScript 'stage_cleanup.sh'

    tag "${sample_id}"

    input:
    tuple val(sample_id), path("tmp/input.bam")
    path("tmp/ref_genome.fasta")
    path("tmp/ref_genome.fasta.fai")
    path("tmp/coverage.bed")

    output:
    tuple val(sample_id), path("vardict.bcf"), emit: calling_result

    script:
    """
    mkdir -p tmp/bed_chunks/
    split -n l/${task.cpus} --additional-suffix=".bed" -a 4 -d tmp/coverage.bed tmp/bed_chunks/

    samtools index --threads ${task.cpus} --csi tmp/input.bam

    mkdir -p tmp/calling_chunks/

    cat <<EOF > tmp/calling_chunks/sample_id.txt
    ${sample_id}
    EOF

    parallel -j ${task.cpus} \
    'vardict-java ${params.calling.vardict.extra_args} -G tmp/ref_genome.fasta -b tmp/input.bam -fisher -th 1 -q "${params.calling.min_base_quality}" -VS SILENT -c 1 -S 2 -E 3 -g 4 {} | var2vcf_valid.pl -S -q "${params.calling.min_base_quality}" -E | bcftools reheader -s tmp/calling_chunks/sample_id.txt -f tmp/ref_genome.fasta.fai | bcftools filter -e"QUAL<${params.calling.min_snp_qual}" -Ou | bcftools norm --check-ref x -m+ -Ou --fasta-ref tmp/ref_genome.fasta -o tmp/calling_chunks/{#}.bcf
    bcftools index --threads 1 --csi tmp/calling_chunks/{#}.bcf' ::: tmp/bed_chunks/*.bed

    bcftools concat --allow-overlaps --threads ${task.cpus} -Ou tmp/calling_chunks/*.bcf | bcftools sort -Ou -o vardict.bcf
    """
}