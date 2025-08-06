#!/bin/bash
set -euo pipefail

echo "Prepearing genome index..."
bowtie2-build --quiet --threads 4 test_data/arth_chr1.fasta test_data/arth_chr1.fasta
echo "Done"

name="sample1"
read1=$(readlink -f test_data/SRR33243472_1.fastq.gz)
read2=$(readlink -f test_data/SRR33243472_2.fastq.gz)
echo -e "${name}\t${read1}\t${read2}\t${read1}" > test_data/samples.tsv
for i in 'sample2' 'sample3'; do echo -e "${i}\t${read1}\t${read2}\t${read1}" >> test_data/samples.tsv; done

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
