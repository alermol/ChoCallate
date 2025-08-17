#!/bin/bash
set -euo pipefail

echo "Prepearing genome index..."
bowtie2-build --quiet --threads 4 test_data/arth_chr1.fasta.gz test_data/arth_chr1.fasta.gz
echo "Done"

name="sample1"
read1=$(readlink -f test_data/test_reads_R1.fq.gz)
read2=$(readlink -f test_data/test_reads_R2.fq.gz)
read3=$(readlink -f test_data/test_reads_SE.fq.gz)
echo -e "${name}\t${read1}\t${read2}\t${read3}" > test_data/samples.tsv

nextflow run main.nf \
    --samples_tsv test_data/samples.tsv \
    --reference_index test_data/arth_chr1.fasta.gz \
    --outdir test_data/chocallate_test \
    --reference_genome test_data/arth_chr1.fasta.gz

