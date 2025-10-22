#!/bin/bash
set -euo pipefail

echo "Preparing genome index..."
if [[ -e test_data/arth_chr1.fasta.gz.1.bt2 ]]; then
  echo "Bowtie2 index already exists, skipping index building."
else
  bowtie2-build --quiet --threads 4 test_data/arth_chr1.fasta.gz test_data/arth_chr1.fasta.gz
  echo "Bowtie2 index built."
fi
echo "Done"

name="sample1"
read1=$(readlink -f test_data/test_reads_R1.fq.gz)
read2=$(readlink -f test_data/test_reads_R2.fq.gz)
read3=$(readlink -f test_data/test_reads_SE.fq.gz)
bam_path=$(readlink -f test_data/sample1.bam)

# TEST 1: input_format=fastq, reads_type=pe (should PASS)
# samples.tsv: name, R1, R2, '-'
echo -e "${name}\t${read1}\t${read2}\t-" > test_data/samples.tsv
echo "==== TEST 1: input_format=fastq, reads_type=pe (should PASS) ===="
nextflow run main.nf \
    --samples_tsv test_data/samples.tsv \
    --input_format fastq \
    --reference_index test_data/arth_chr1.fasta.gz \
    --outdir test_data/chocallate_test_1 \
    --reference_genome test_data/arth_chr1.fasta.gz || echo "FAILED as expected?"

# TEST 2: input_format=bam, reads_type=pe (should PASS)
# samples.tsv: name, bam_path
echo -e "${name}\t${bam_path}" > test_data/samples.tsv
echo "==== TEST 2: input_format=bam, reads_type=pe (should PASS) ===="
nextflow run main.nf \
    --samples_tsv test_data/samples.tsv \
    --input_format bam \
    --outdir test_data/chocallate_test_2 \
    --reference_genome test_data/arth_chr1.fasta.gz || echo "FAILED as expected?"

