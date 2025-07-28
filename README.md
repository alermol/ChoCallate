# ChoCallate 🍫

**ChoCallate** (**Cho**rus of **Call**ers) is an automated pipeline for calling single-nucleotide variants (SNVs) and indels (INDELs) using several variant callers. The pipeline consolidates results and applies a majority rule for final variant calling.

## 📋 Description
<!-- 
⚠️ **Important Notes**
- **Diploid species**: The main pipeline `ChoCallate.nf` is optimized for diploid organisms. Some integrated callers (e.g., vardict, bcftools) assume diploidy by default
- **Polyploid species**: For polyploid species, use the `polyChoCallate.nf` pipeline with explicit ploidy parameterization. VarDict and bcftools are not using in `polyChoCallate.nf`. -->


**Key Features**
- **Consensus-driven approach**: Combines results from 5 callers for diploid calling (FreeBayes, bcftools, GATK4, VarDict, SNVer) or 3 for polyploid calling (VarDict and bcftools are not support ploidy customization), using majority rule
- **Ploidy flexibility**: Supports both diploid and polyploid species.
- **Flexible input**: Compatible with both GBS (Genotyping-by-Sequencing) and WGS data
- **Quality filtering**: Multiple filtering steps based on coverage, base quality, and SNP quality


**Workflow Steps**
1. **Alignment**: Bowtie2-based read alignment
2. **Variant Calling**: Parallel execution of 5 or 3 callers
3. **Consensus Generation**: Merges results using majority rule
4. **Output**: Final compressed VCF files

![image](ChoCallate_scheme.png)


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
bash run_test.sh
```
4. **(Optional) Cleanup output after test**
```bash
bash cleanup.sh
```
5. **(Optional) Adjust amout of RAM available for Snver and VarDict**     
Important if you are working with large genomes or large numbers of reads      
Replace N with the desired amount of RAM in GB
```bash
sed -i 's/-Xmx1g/-XmxNg/' $CONDA_PREFIX/bin/snver
sed -i 's/-Xmx8g/-XmxNg/' $CONDA_PREFIX/bin/vardict-java
```      


## 🚀 Usage
**Basic execution**:
```bash
./ChoCallate.py \
    --reference_genome /path/to/ref.fasta \
    --reference_index /path/to/ref_index \
    --samples_tsv samples.tsv
```

## ⚙️ Parameters
| Parameter | Default | Description | Didploid calling | Polyploid calling |
| :-------- | :------- | :---------- | :---------: | :----------: |
|--samples_tsv	|samples.tsv|	TSV file with samples| 🟢 | 🟢 |
|--outdir	|ChoCallate_output|	Output directory| 🟢 | 🟢 |
|--reference_genome	|-|	Reference genome (FASTA), must be not gzipped| 🟢 | 🟢 |
|--reference_index	|-	|Bowtie2 index prefix for the reference| 🟢 | 🟢 |
|--min_coverage|	5|	Minimum position coverage depth for SNP-calling| 🟢 | 🟢 |
|--min_base_quality	|20	|Minimum base quality for SNP-calling| 🟢 | 🟢 |
|--min_snp_qual	|20	|Minimum SNP quality| 🟢 | 🟢 |
|--reads_type	|pe	|Reads type (pe for paired-end, se for single-end)| 🟢 | 🟢 |
|--reads_source	|gbs|	Data source (gbs for GBS or wgs for WGS)| 🟢 | 🟢 |
|--bowtie2_cpu	|10	|Number of threads for Bowtie2| 🟢 | 🟢 |
|--bowtie2_forks	|1	|Number of Bowtie2 processes running in parallel| 🟢 | 🟢 |
|--freebayes_forks	|1	|Number of freebayes processes running in parallel| 🟢 | 🟢 |
|--gatk4_forks	|1	|Number of GATK4 processes running in parallel| 🟢 | 🟢 |
|--snver_forks	|1	|Number of Snver processes running in parallel| 🟢 | 🟢 |
|--snps_cons_forks	|1	|Number of SNPs consensus VCF generation processes running in parallel| 🟢 | 🟢 |
|--indels_cons_forks	|1	|Number of INDELs consensus VCF sgeneration processes running in parallel| 🟢 | 🟢 |
|--chunk_size	|0	|Chunk size for splitting the input file list. A value of 0 means file will not be split into chunks. For details, see the footnote below the table.| 🟢 | 🟢 |
|--debug	|false	|If set, the working directory will not be deleted after the pipeline completes| 🟢 | 🟢 |
|--bcftools_cpu	|1	|Number of threads for bcftools| 🟢 | 🔴 |
|--bcftools_forks	|1	|Number of bcftools processes running in parallel| 🟢 | 🔴 |
|--vardict_cpu	|1	|Number of threads for VarDict| 🟢 | 🔴 |
|--vardict_forks	|1	|Number of VarDict processes running in parallel| 🟢 | 🔴 |
|--ploidy |2 |Ploidy | 🔴 | 🟢 |

> Each chunk will be processed sequentially. Intermediate files will be deleted after processing each chunk is complete. This is useful when processing a large number of files, which could otherwise lead to running out of disk space due to temporary files. 


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
Ermolaev, A. (2025). *ChoCallate: Consensus variant calling pipeline* \[Computer software\]. GitHub. https://github.com/alermol/ChoCallate

**BibTeX**:  
```bibtex
@software{ChoCallate,
  author = {Ermolaev, A.},
  title = {ChoCallate: Consensus variant calling pipeline},
  url = {https://github.com/alermol/ChoCallate},
  year = {2025}
}
```

## 📜 License
MIT License
