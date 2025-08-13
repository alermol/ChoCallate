# ChoCallate 🍫

**ChoCallate** (**Cho**rus of **Call**ers) is an automated pipeline for calling single-nucleotide variants (SNVs) and indels (INDELs) using several variant callers. The pipeline consolidates results and applies configurable consensus rules for final variant calling.

## 📋 Description

**Key Features**
- **Consensus-driven approach**: Combines results from multiple variant callers using configurable consensus rules
- **Multiple consensus types**: Supports majority rule (mj), n-1 consensus (n1), and full consensus (fc)
- **Ploidy flexibility**: Supports both diploid and polyploid species with automatic caller selection
- **Flexible input**: Compatible with both GBS (Genotyping-by-Sequencing) and WGS data
- **Quality filtering**: Multiple filtering steps based on coverage, base quality, and SNP quality
- **Parallel processing**: Efficient parallel execution of variant callers and consensus generation

**Supported Variant Callers**
- **Diploid calling**: bcftools, GATK4, FreeBayes, SNVer, VarDict
- **Polyploid calling**: GATK4, FreeBayes, SNVer (bcftools and VarDict don't support ploidy customization)

**Workflow Steps**
1. **Alignment**: Bowtie2-based read alignment with quality filtering
2. **Variant Calling**: Parallel execution of selected variant callers
3. **Consensus Generation**: Merges results using configurable consensus rules
4. **Output**: Final compressed VCF files for SNPs and INDELs

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

3. **Run the pipeline on test data**:
```bash
bash run_test.sh
```

4. **(Optional) Cleanup output after test**:
```bash
bash cleanup.sh
```

5. **(Optional) Adjust amount of RAM available for SNVer and VarDict**     
Important if you are working with large genomes or large numbers of reads      
Replace N with the desired amount of RAM in GB:
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

**Advanced execution with custom consensus and ploidy**:
```bash
./ChoCallate.py \
    --reference_genome /path/to/ref.fasta \
    --reference_index /path/to/ref_index \
    --samples_tsv samples.tsv \
    --ploidy 4 \
    --cons_type n1 \
    --effective_callers gatk,freebayes,snver
```

## ⚙️ Parameters

### Required Parameters
| Parameter | Description |
| :-------- | :---------- |
| `--reference_genome` | Reference genome in FASTA format (must not be gzipped) |
| `--reference_index` | Bowtie2 index prefix for the reference genome |

### Input/Output Parameters
| Parameter | Default | Description |
| :-------- | :------- | :---------- |
| `--samples_tsv` | `samples.tsv` | TSV file with sample information |
| `--outdir` | `ChoCallate_output` | Output directory for results |
| `--chunk_size` | `0` | Chunk size for splitting input files (0 = no splitting) |

### Quality and Filtering Parameters
| Parameter | Default | Description |
| :-------- | :------- | :---------- |
| `--min_coverage` | `5` | Minimum position coverage depth for variant calling |
| `--min_base_quality` | `20` | Minimum base quality for variant calling |
| `--min_map_qual` | `10` | Minimum mapping quality for read filtering |
| `--min_snp_qual` | `20` | Minimum variant quality threshold |

### Data Type Parameters
| Parameter | Default | Choices | Description |
| :-------- | :------- | :------ | :---------- |
| `--reads_type` | `pe` | `pe`, `se`, `mx` | Read type: paired-end, single-end, or mixed |
| `--reads_source` | `gbs` | `gbs`, `wgs` | Data source: GBS or whole genome sequencing |
| `--ploidy` | `2` | `≥2` | Ploidy level of the organism |

### Variant Calling Parameters
| Parameter | Default | Description |
| :-------- | :------- | :---------- |
| `--effective_callers` | `-` | Comma-separated list of variant callers to use. Use `-` for automatic selection based on ploidy |
| `--cons_type` | `mj` | Consensus type: `mj` (majority), `n1` (n-1), `fc` (full consensus) |

### Resource Allocation Parameters
| Parameter | Default | Description |
| :-------- | :------- | :---------- |
| `--bowtie2_cpu` | `10` | Number of threads for Bowtie2 alignment |
| `--bowtie2_forks` | `1` | Number of parallel Bowtie2 processes |
| `--zero_vcf_cpu` | `1` | Number of threads for zero VCF generation |
| `--zero_vcf_forks` | `1` | Number of parallel zero VCF processes |
| `--cons_cpus` | `5` | Number of threads for consensus generation |
| `--cons_forks` | `1` | Number of parallel consensus processes |
| `--bcftools_cpu` | `1` | Number of threads for bcftools |
| `--vardict_cpu` | `1` | Number of threads for VarDict |

### Processing Parameters
| Parameter | Default | Description |
| :-------- | :------- | :---------- |
| `--win_size` | `1000000` | Window size (in bp) for parallel consensus generation |
| `--debug` | `false` | Keep working directory after pipeline completion |

### Consensus Types
- **`mj` (Majority Rule)**: Variant is called if majority of callers identify it
- **`n1` (N-1 Consensus)**: Variant is called if n-1 callers identify it (where n is total number of callers)
- **`fc` (Full Consensus)**: Variant is called only if all callers identify it

> **Note**: Each chunk will be processed sequentially. Intermediate files will be deleted after processing each chunk is complete. This is useful when processing a large number of files, which could otherwise lead to running out of disk space due to temporary files.

### Automatic Caller Selection
When `--effective_callers` is set to `-` (default), ChoCallate automatically selects appropriate callers among available:
- **Diploid (ploidy=2)**: Uses `bcftools,gatk,freebayes,snver,vardict`
- **Polyploid (ploidy>2)**: Uses `gatk,freebayes,snver` (polyploid-compatible callers only)

## 📂 Input Data Structure

### Samples TSV Format
Example `samples.tsv`:
```text
sample1    /path/to/sample1_R1.fq.gz    /path/to/sample1_R2.fq.gz    /path/to/sample1_SE.fq.gz
sample2    /path/to/sample1_R1.fq.gz    /path/to/sample1_R2.fq.gz    /path/to/sample1_SE.fq.gz
```

**`samples.tsv` requirements**:
- Field separator must be tab
- All paths must be absolute
- Format: `sample_id`, `read1`, `read2`, `read3`

**Read type configurations**:
- **`--reads_type pe`**: Use columns 2 and 3 for paired-end reads, column 4 can be any symbol
- **`--reads_type se`**: Use column 2 for single-end reads, columns 3 and 4 should contain the same data
- **`--reads_type mx`**: Use columns 2 and 3 for paired-end reads, column 4 for single-end reads

**Reference requirements**:
- Reference genome in FASTA format (uncompressed)
- Bowtie2 index must be pre-built

⚠️ **Important Note**     
Since ChoCallate uses Bowtie2 for alignment, it is possible to process multiple archives of reads simultaneously using comma-separated lists. See the [Bowtie2 documentation](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#paired-inputs:~:text=BOWTIE2_INDEXES%20environment%20variable.-,%2D1%20%3Cm1%3E,bowtie2%20gets%20the%20reads%20from%20the%20%22standard%20in%22%20or%20%22stdin%22%20filehandle.,-%2D%2Dinterleaved) for details.

## 📊 Output

The `outdir` will contain:
- `<sample>` folder for each sample
- `<sample>/<sample>.snps.vcf.gz` — Final VCF containing SNPs only (bgzipped)
- `<sample>/<sample>.indels.vcf.gz` — Final VCF containing INDELs only (bgzipped)
- Pipeline reports (HTML format) and execution traces

## 🛠️ Dependencies

All dependencies are managed via Conda (`environment.yaml`):
- **Variant Callers**: FreeBayes, bcftools, GATK4, SNVer, VarDict
- **Alignment**: Bowtie2
- **Processing**: samtools, bedtools, bedops
- **Pipeline**: Nextflow
- **Utilities**: Python, tabix, parallel

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
