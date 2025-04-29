#!/bin/bash
set -euo pipefail

echo "Prepearing genome index..."
bowtie2-build --quiet --threads 4 arth_chr1.fasta arth_chr1.fasta
echo "Done"

name="test"
read1=$(readlink -f ERR1839886_1.fastq.gz)
read2=$(readlink -f ERR1839886_2.fastq.gz)
echo -e "${name}\t${read1}\t${read2}" > samples.tsv

nextflow run ../ChoCallate.nf \
    -c ../ChoCallate.config \
    --samples_tsv samples.tsv \
    --reference_index arth_chr1.fasta \
    --outdir chocallate_test \
    --reference_genome arth_chr1.fasta
