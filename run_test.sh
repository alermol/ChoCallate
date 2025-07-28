#!/bin/bash
set -euo pipefail

echo "Prepearing genome index..."
bowtie2-build --quiet --threads 4 test_data/arth_chr1.fasta test_data/arth_chr1.fasta
echo "Done"

name="test"
read1=$(readlink -f test_data/SRR33243472_1.fastq.gz)
read2=$(readlink -f test_data/SRR33243472_2.fastq.gz)
echo -e "${name}\t${read1}\t${read2}" > test_data/samples.tsv

./ChoCallate.py \
    --samples_tsv test_data/samples.tsv \
    --reference_index test_data/arth_chr1.fasta \
    --outdir test_data/chocallate_test \
    --reference_genome test_data/arth_chr1.fasta

./ChoCallate.py \
    --samples_tsv test_data/samples.tsv \
    --reference_index test_data/arth_chr1.fasta \
    --outdir test_data/chocallate_test \
    --reference_genome test_data/arth_chr1.fasta \
    --ploidy 4
