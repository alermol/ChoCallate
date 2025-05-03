# ChoCallate 🍫

**ChoCallate** (*Chor*us of *Call*ers) is an automated pipeline for calling single-nucleotide variants (SNVs) and indels (INDELs) using multiple popular variant callers. The pipeline consolidates results and applies a majority rule for final variant calling.

## 📋 Description

⚠️ **Important Note**  
This pipeline is designed specifically for diploid organisms. Some integrated callers (e.g., vardict, bcftools) assume diploidy by default and lack explicit ploidy parameterization.

**Key Features**
- **Consensus-driven approach**: Combines results from 5 callers (FreeBayes, bcftools, GATK4, VarDict, SNVer) using majority rule (≥3/5 callers)
- **Flexible input**: Compatible with both GBS (Genotyping-by-Sequencing) and WGS data
- **Quality filtering**: Multiple filtering steps based on coverage, base quality, and SNP quality

**Workflow Steps**
1. **Alignment**: Bowtie2-based read alignment
2. **Variant Calling**: Parallel execution of 5 callers
3. **Consensus Generation**: Merges results using majority rule
4. **Output**: Final compressed VCF files


## ⚙️ Installation

1. **Clone the repository**:
```bash
git clone https://github.com/alermol/ChoCallate.git
cd ChoCallate
```
2. **Set up the Conda environment**:
```bash
conda env create -f environment.yaml
conda activate ChoCallate
```
3. **Run the pipeline on a test data**:
```bash
cd test_data
bash run_test.sh
```
4. **(Optional) Cleanup output after test**
```bash
bash cleanup.sh
```

## 🚀 Usage
**Basic execution**:
```bash
nextflow run ChoCallate.nf \
    -c ChoCallate.config \
    --reference_genome /path/to/ref.fasta \
    --reference_index /path/to/ref_index \
    --samples_tsv samples.tsv
```

## ⚙️ Parameters
| Parameter | Default | Description |
| :-------- | :------- | :---------- |
|--samples_tsv	|samples.tsv|	TSV file with samples|
|--outdir	|ChoCallate_output|	Output directory|
|--reference_genome	|-|	Reference genome (FASTA), must be not gzipped|
|--reference_index	|-	|Bowtie2 index prefix for the reference|
|--min_coverage|	5|	Minimum position coverage depth for SNP-calling|
|--min_base_quality	|20	|Minimum base quality for SNP-calling|
|--bowtie2_cpu	|10	|Number of threads for Bowtie2|
|--min_snp_qual	|20	|Minimum SNP quality|
|--reads_type	|pe	|Reads type (pe for paired-end, se for single-end)|
|--reads_source	|gbs|	Data source (gbs for GBS or wgs for WGS)|


## 📂 Input Data Structure
Example `samples.tsv`:    
```text
sample1    /path/to/sample1_R1.fq.gz    /path/to/sample1_R2.fq.gz
sample2    /path/to/sample2_R1.fq.gz    /path/to/sample2_R2.fq.gz
````

**`samples.tsv` requirements**:
- Field separator for `samples.tsv` must be tab
- All paths in `samples.tsv` must be absolute

**Reference requirements**:
- Reference genome in FASTA format
- Bowtie2 index must be pre-built

If `--reads_type se`, then the third column in `samples.tsv` must contain the same data as the second one.


## 📊 Output
The `outdir` will contain:
- `<sample>` folder
- `<sample>`/`<sample>`.snps.vcf.gz — Final VCF for the `<sample>` containing SNPs only (bgziped)
- `<sample>`/`<sample>`.indels.vcf.gz — Final VCF for the `<sample>` containing INDELs only (bgziped)


## 🛠️ Dependencies
All dependencies are managed via Conda (`environment.yaml`)


## ❓ Support
For questions or issues, please open an issue on GitHub.


## ✏️ Citation

To cite this software repository in your work, use the following format:

**APA Style**:  
Ermolaev, A. (2025). *ChoCallate: Consensus variant calling pipeline for diploid organisms* \[Computer software\]. GitHub. https://github.com/alermol/ChoCallate

**BibTeX**:  
```bibtex
@software{ChoCallate,
  author = {Ermolaev, A.},
  title = {ChoCallate: Consensus variant calling pipeline for diploid organisms},
  url = {https://github.com/alermol/ChoCallate},
  year = {2025}
}
```

## 📜 License
MIT License
