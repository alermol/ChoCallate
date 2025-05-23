// ChoCallate
// Pipeline for SNV/INDEL calling using *Cho*rus of *Call*ers

params.samples_tsv = 'input.tsv'
params.outdir = 'ChoCallate_output'
params.min_coverage = 5
params.min_base_quality = 20
params.bowtie2_cpu = 10
params.samtools_min_map_qual = 10
params.min_snp_qual = 20
params.reads_type = 'pe' // se - single-end reads; pe - pair-end reads
params.reads_source = 'gbs' // gbs - Genotyping-by-sequencing; wgs - Whole Genome Sequencing


// Main workflow definition
workflow {
    // Create channel from samples TSV file and parse into tuples
    Channel
        .fromPath( params.samples_tsv )
        .splitCsv( header: false, sep: '\t' )
        .map { row -> tuple( row[0], file(row[1]), file(row[2]) ) }
        .set { sample_run_ch }

    // Define reference files
    ref_index = file(params.reference_index)
    ref_genome = file(params.reference_genome)

    // Create reference genome index files
    create_faidx(ref_genome)
    create_sequence_dictionary(ref_genome)

    // Map reads to reference genome
    map_reads(sample_run_ch, ref_index)

    // Left-align indels in BAM files
    left_align_indels(map_reads.out.bam, 
                      create_sequence_dictionary.out.gen_dict, 
                      create_faidx.out.ref_genome)

    // Index the aligned BAM files
    bam_indexing(left_align_indels.out.lai_bam)

    // Generate coverage information from BAM files
    bam_cov_generation(bam_indexing.out.ind_bam.map{it[0]})

    // Call SNVs/INDELs using different callers
    freebayes_snps = freebayes_calling(bam_indexing.out.ind_bam, 
                                      bam_cov_generation.out.coverage, 
                                      create_faidx.out.ref_genome)
    bcftools_snps = bcftools_calling(bam_indexing.out.ind_bam, 
                                    bam_cov_generation.out.coverage, 
                                    create_faidx.out.ref_genome)
    gatk4_snps = gatk4_calling(bam_indexing.out.ind_bam, 
                              bam_cov_generation.out.coverage, 
                              create_faidx.out.ref_genome,
                              create_sequence_dictionary.out.gen_dict)
    vardict_snps = vardict_calling(bam_indexing.out.ind_bam, 
                                  bam_cov_generation.out.coverage, 
                                  create_faidx.out.ref_genome)
    snver_snps = snver_calling(bam_indexing.out.ind_bam, 
                              bam_cov_generation.out.coverage, 
                              create_faidx.out.ref_genome)

    // Process SNPs
    merged_snps_vcfs = freebayes_snps.snps_vcf
        .join(bcftools_snps.snps_vcf)
        .join(gatk4_snps.snps_vcf)
        .join(vardict_snps.snps_vcf)
        .join(snver_snps.snps_vcf)
        .map { sample, freebayes, bcftools, gatk, vardict, snver -> 
            tuple(sample, freebayes, bcftools, gatk, vardict, snver) 
        }
    generate_final_vcf_snps(merged_snps_vcfs)
    process_final_vcf_snps(generate_final_vcf_snps.out.fvcf, create_faidx.out.ref_genome.map{it[1]})

    // Process INDELs
    merged_indels_vcfs = freebayes_snps.indels_vcf
        .join(bcftools_snps.indels_vcf)
        .join(gatk4_snps.indels_vcf)
        .join(vardict_snps.indels_vcf)
        .join(snver_snps.indels_vcf)
        .map { sample, freebayes, bcftools, gatk, vardict, snver -> 
            tuple(sample, freebayes, bcftools, gatk, vardict, snver) 
        }
    generate_final_vcf_indels(merged_indels_vcfs)
    process_final_vcf_indels(generate_final_vcf_indels.out.fvcf, create_faidx.out.ref_genome.map{it[1]})
}

// Cleanup temporary files after workflow completion
workflow.onComplete {
    def tmpDir = file("work/")
    if (tmpDir.exists()) {
        tmpDir.deleteDir()
        }
}

// Process to create FASTA index file (for FreeBayes)
process create_faidx {
    maxForks 1
    cpus 1
    
    tag "${ref_genome}-faidx"

    input:
    path(ref_genome)

    output:
    tuple path("${ref_genome}"), path("${ref_genome}.fai"), emit: ref_genome

    script:
    """
    samtools faidx --threads ${task.cpus} ${ref_genome}
    """
}

// Process to create sequence dictionary (for GATK4)
process create_sequence_dictionary {
    tag "${ref_genome.baseName}-CreateSequenceDictionary"

    input:
    path(ref_genome)
    
    output:
    path("${ref_genome.baseName}.dict"), emit: gen_dict
    
    script:
    """
    gatk CreateSequenceDictionary -R ${ref_genome}
    """
}

// Process to map reads to reference genome
process map_reads {
    maxForks 1
    cpus params.bowtie2_cpu

    tag "${sample_id}-bowtie2"
    
    input:
    tuple val(sample_id), val(read1), val(read2)
    val(ref_index)

    output:
    path("${sample_id}.tmp_bam"), emit: bam

    script:
    if ( params.reads_type == 'pe' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -1 ${read1} -2 ${read2} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.samtools_min_map_qual} -F 4 - | \
            samtools fixmate -@ ${task.cpus} -m - - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.tmp_bam
        """
    else if ( params.reads_type == 'se' )
        """
        bowtie2 --threads ${task.cpus} --rg-id ${sample_id} --rg SM:${sample_id} -x ${ref_index} -U ${read1} | \
            samtools view -@ ${task.cpus} -S -b -q ${params.samtools_min_map_qual} -F 4 - | \
            samtools sort -@ ${task.cpus} -o ${sample_id}.tmp_bam
        """
    else
        error 'Invalid reads type: ${params.reads_type}. Available types: se, pe'
}

// Process to left-align indels in BAM files
process left_align_indels {
    cpus 1
    maxForks 1

    tag "${bam.baseName}-LeftAlignIndels"

    input:
    path(bam)
    path(genome_dictionary)
    tuple path(ref_genome), path(genome_fai)
    
    output:
    path("${bam.baseName}.bam"), emit: lai_bam
    
    script:
    """
    gatk LeftAlignIndels -I ${bam} -O ${bam.baseName}.bam -R ${ref_genome} -OBI false
    """
}

// Process to index BAM files
process bam_indexing {
    maxForks 1
    cpus 1

    tag "${bam.baseName}-samtools"
    
    input:
    path(bam)
    
    output:
    tuple path("${bam}"), path("${bam}.csi"), emit: ind_bam
    
    script:
    """
    samtools index --csi --threads ${task.cpus} ${bam.baseName}.bam
    """
}

// Process to generate coverage information
process bam_cov_generation {
    cpus 1
    maxForks 1
    
    tag "${bam.baseName}-coverage"

    input:
    path(bam)

    output:
    path("${bam.baseName}.bed"), emit: coverage

    script:
    """
    samtools depth --threads ${task.cpus} ${bam} | \
        awk '\$3 >= ${params.min_coverage} {print \$1,\$2-1,\$2}' | \
        bedops --merge - > ${bam.baseName}.bed
    """
}

// Process to call variants using FreeBayes
process freebayes_calling {
    cpus 1

    tag "${bam.baseName}-freebayes"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_freebayes"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_freebayes"), emit: indels_vcf
    
    script:
    if ( params.reads_source == 'gbs' )
        """
        freebayes --fasta-reference ${ref_genome} --targets ${coverage} --dont-left-align-indels \
            --use-best-n-alleles 4 --min-alternate-qsum ${params.min_base_quality} --hwe-priors-off --no-population-priors \
            --binomial-obs-priors-off --allele-balance-priors-off --min-base-quality ${params.min_base_quality} \
            --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam ${bam} --limit-coverage 250 | \
            bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
            bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
            bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.freebayes

        bcftools view -v snps -Ov -o ${bam.baseName}.snps_freebayes ${bam.baseName}.freebayes
        
        bcftools view -v indels -Ov -o ${bam.baseName}.indels_freebayes ${bam.baseName}.freebayes
        """
    
    else if ( params.reads_source == 'wgs' )
        """
        freebayes --fasta-reference ${ref_genome} --targets ${coverage} --dont-left-align-indels \
            --use-best-n-alleles 4 --min-alternate-qsum ${params.min_base_quality} --hwe-priors-off --no-population-priors \
            --allele-balance-priors-off --min-base-quality ${params.min_base_quality} \
            --haplotype-length -1 --throw-away-complex-obs --no-partial-observations --bam ${bam} --limit-coverage 250 | \
            bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
            bcftools view --min-alleles 2 --max-alleles 2 - | bcftools annotate --force -x INFO,FORMAT - | \
            bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.freebayes

        bcftools view -v snps -Ov -o ${bam.baseName}.snps_freebayes ${bam.baseName}.freebayes
        
        bcftools view -v indels -Ov -o ${bam.baseName}.indels_freebayes ${bam.baseName}.freebayes
        """
    
    else
        error 'Invalid reads source: ${params.reads_source}. Available sources: gbs, wgs'
}

// Process to call variants using bcftools
process bcftools_calling {
    cpus 1
    
    tag "${bam.baseName}-bcftools"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_bcftools"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_bcftools"), emit: indels_vcf
    
    script:
    """
    bcftools mpileup -Ou --count-orphans --fasta-ref ${ref_genome} --threads ${task.cpus} --max-depth 250 \
        --min-BQ ${params.min_base_quality} --regions-file ${coverage} ${bam} | \
        bcftools call -Ov --multiallelic-caller --threads ${task.cpus} | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.bcftools

    bcftools view -V indels,mnps,bnd,other -Ov -o ${bam.baseName}.snps_bcftools ${bam.baseName}.bcftools
    
    bcftools view -v indels -Ov -o ${bam.baseName}.indels_bcftools ${bam.baseName}.bcftools
    """
}

// Process to call variants using GATK4
process gatk4_calling {
    cpus 1

    tag "${bam.baseName}-gatk4"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    path(ref_genome_dict)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_gatk"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_gatk"), emit: indels_vcf
    
    script:
    """
    gatk HaplotypeCaller -I ${bam} -R ${ref_genome} -mbq ${params.min_base_quality} -O ${bam.baseName}.gatk1 -L ${coverage}
    bcftools filter ${bam.baseName}.gatk1 -e'QUAL<${params.min_snp_qual}' | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools sort - | \
        bcftools view -AA --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.gatk

    bcftools view -v snps -Ov -o ${bam.baseName}.snps_gatk ${bam.baseName}.gatk
    
    bcftools view -v indels -Ov -o ${bam.baseName}.indels_gatk ${bam.baseName}.gatk
    """
}

// Process to call variants using VarDict
process vardict_calling {
    cpus 1

    tag "${bam.baseName}-vardict"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_vardict"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_vardict"), emit: indels_vcf
    
    script:
    """
    vardict-java -G ${ref_genome} -N ${bam.baseName} -b ${bam} -fisher -th ${task.cpus} \
        -VS SILENT --nosv -k 0 -q ${params.min_base_quality} -c 1 -S 2 -E 3 -g 4 ${coverage} | \
        var2vcf_valid.pl -q ${params.min_base_quality} -N ${bam.baseName} -E | \
        bcftools reheader -f ${ref_genome_fai} - | bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.vardict

    bcftools view -v snps -Ov -o ${bam.baseName}.snps_vardict ${bam.baseName}.vardict
    
    bcftools view -v indels -Ov -o ${bam.baseName}.indels_vardict ${bam.baseName}.vardict
    """
}

// Process to call variants using SNVer
process snver_calling {
    cpus 1

    tag "${bam.baseName}-snver"
    
    input:
    tuple path(bam), path(bam_index)
    path(coverage)
    tuple path(ref_genome), path(ref_genome_fai)
    
    output:
    tuple val("${bam.baseName}"), path("${bam.baseName}.snps_snver"), emit: snps_vcf
    tuple val("${bam.baseName}"), path("${bam.baseName}.indels_snver"), emit: indels_vcf
    
    script:
    """
    ln -s ${ref_genome} reference.fasta
    
    samtools faidx reference.fasta
    
    snver -i ${bam} -r reference.fasta -o ${bam.baseName} -l ${coverage} -bq ${params.min_base_quality}
    
    bcftools reheader -f reference.fasta.fai ${bam.baseName}.filter.vcf | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.snps_snver
        
    bcftools reheader -f reference.fasta.fai ${bam.baseName}.indel.filter.vcf | \
        bcftools filter -e'QUAL<${params.min_snp_qual}' - | \
        bcftools annotate --force -x INFO,FORMAT - | bcftools view --min-alleles 2 --max-alleles 2 - | \
        bcftools norm --fasta-ref ${ref_genome} --atom-overlaps '.' --atomize -Ov -o ${bam.baseName}.indels_snver
    """
}


// Process to generate final consensus VCF from all callers using majority rule
// Majority rule - variant is true if detected more than 2 callers
process generate_final_vcf_indels {
    maxForks 1
    cpus 1

    tag "${sample}-generate"

    input:
    tuple val(sample), path(vcf1), path(vcf2), path(vcf3), path(vcf4), path(vcf5)

    output:
    path("${sample}.vcf"), emit: fvcf

    script:
    """
    #!/usr/bin/env python3
    
    import sys
    from collections import defaultdict, Counter

    def sort_gt(s):
        parts = s.split('/')
        numbers = [int(part) for part in parts]
        numbers_sorted = sorted(numbers)
        result = '/'.join(map(str, numbers_sorted))
        return result

    
    def get_most_frequent(items):
        counter = Counter(items)
        most_common = counter.most_common(1)
        return most_common
    
    
    def parse_vcf(vcf_file):
        variants = defaultdict(dict)
        with open(vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\\t')
                if '.' in fields[9]:
                    continue
                else:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    gt = sort_gt(fields[9])
                    variants[(chrom, pos)] = [(ref, alt, gt)]
        return variants

    indels = defaultdict(list)
    for file in ["${vcf1}", "${vcf2}", "${vcf3}", "${vcf4}", "${vcf5}"]:
        polymorphic_indels = parse_vcf(file)
        for coord, gen in polymorphic_indels.items():
            indels[coord].append(gen[0])

    with open("${sample}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\\n')
        out.write(f'#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t{"${sample}"}\\n')
        for coord, gen in indels.items():
            most_freq_var = get_most_frequent(gen)
            if ((most_freq_var[0][1] < 3) | 
                ('.' in most_freq_var[0][0][2]) | 
                (not most_freq_var[0][0][0].isupper())):
                continue
            else:
                chrom = coord[0]
                pos = coord[1]
                ref = most_freq_var[0][0][0]
                alt = most_freq_var[0][0][1]
                gt = most_freq_var[0][0][2]
                out.write('\\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', '.', 'GT', f'{gt}\\n']))
    """
}

process generate_final_vcf_snps {
    maxForks 1
    cpus 1

    tag "${sample}-generate"

    input:
    tuple val(sample), path(vcf1), path(vcf2), path(vcf3), path(vcf4), path(vcf5)

    output:
    path("${sample}.vcf"), emit: fvcf

    script:
    """
    #!/usr/bin/env python3
    
    import sys
    from collections import defaultdict, Counter

    def sort_gt(s):
        parts = s.split('/')
        numbers = [int(part) for part in parts]
        numbers_sorted = sorted(numbers)
        result = '/'.join(map(str, numbers_sorted))
        return result

    
    def get_most_frequent(items):
        counter = Counter(items)
        most_common = counter.most_common(1)
        return most_common
    
    
    def parse_vcf(vcf_file):
        variants = defaultdict(dict)
        with open(vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\\t')
                if '.' in fields[9]:
                    continue
                else:
                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]
                    gt = sort_gt(fields[9])
                    variants[(chrom, pos)] = [(ref, alt, gt)]
        return variants


    base_vars = parse_vcf("${vcf2}")

    for file in ["${vcf1}", "${vcf3}", "${vcf4}", "${vcf5}"]:
        polymorphic_vcf = parse_vcf(file)
        for coord, gen in base_vars.items():
            if coord in polymorphic_vcf.keys():
                base_vars[coord].append(polymorphic_vcf[coord][0])
            else:
                base_vars[coord].append(base_vars[coord][0])

    with open("${sample}.vcf", 'w') as out:
        out.write('##fileformat=VCFv4.3\\n')
        out.write('##FORMAT=<ID=GT,Number=1,Type=String>\\n')
        out.write(f'#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\t{"${sample}"}\\n')
        for coord, gen in base_vars.items():
            most_freq_var = get_most_frequent(gen)
            if ((most_freq_var[0][1] < 3) | 
                ('.' in most_freq_var[0][0][2]) | 
                (not most_freq_var[0][0][0].isupper())):
                continue
            else:
                chrom = coord[0]
                pos = coord[1]
                ref = most_freq_var[0][0][0]
                alt = most_freq_var[0][0][1]
                gt = most_freq_var[0][0][2]
                out.write('\\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', '.', 'GT', f'{gt}\\n']))
    """
}

// Process to finalize and compress the VCF file
process process_final_vcf_snps {
    maxForks 1
    cpus 1

    tag "${vcf.baseName}-finalizeSNP"
    publishDir "${params.outdir}/${vcf.baseName}/", mode: 'copy', pattern: '*.snps.vcf.gz'

    input:
    path(vcf)
    path(ref_genome_fai)

    output:
    path("${vcf.baseName}.snps.vcf.gz")

    script:
    """
    bcftools reheader -f ${ref_genome_fai} ${vcf} | bcftools sort -Oz -o ${vcf.baseName}.snps.vcf.gz
    """
}

process process_final_vcf_indels {
    maxForks 1
    cpus 1

    tag "${vcf.baseName}-finalizeINDEL"
    publishDir "${params.outdir}/${vcf.baseName}/", mode: 'copy', pattern: '*.indels.vcf.gz'

    input:
    path(vcf)
    path(ref_genome_fai)

    output:
    path("${vcf.baseName}.indels.vcf.gz")

    script:
    """
    bcftools reheader -f ${ref_genome_fai} ${vcf} | bcftools sort -Oz -o ${vcf.baseName}.indels.vcf.gz
    """
}
