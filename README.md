# ChoCallate 🍫

**ChoCallate** (**Cho**rus of **Call**ers) is a high-performance, automated pipeline for consensus-based variant calling that combines multiple variant callers to produce robust, high-confidence single-nucleotide variants (SNVs) and indels (INDELs).

## 🎯 What is ChoCallate?

ChoCallate addresses a critical challenge in variant calling: individual variant callers can produce different results for the same genomic data, leading to uncertainty in variant identification. By implementing a **consensus-driven approach**, ChoCallate combines results from multiple state-of-the-art variant callers and applies configurable consensus rules to generate reliable, high-quality variant calls.

### 🌟 Key Features

- **🔄 Consensus-driven approach**: Combines multiple variant callers using configurable consensus rules
- **🧬 Ploidy flexibility**: Supports both diploid and polyploid species with automatic caller selection
- **📊 Multiple consensus types**: Majority rule, n-1 consensus, and full consensus options
- **📂 Dual input support**: Processes both FASTQ (raw reads) and BAM (pre-aligned) files, allowing flexible integration of sequencing data at different analysis stages
- **🔧 Flexible input compatibility**: Works with GBS (Genotyping-by-Sequencing) and WGS data
- **⚡ Parallel processing**: Efficient parallel execution for optimal performance
- **🎛️ Configurable quality filtering**: Multiple filtering steps based on coverage, base quality, and SNP quality
- **📈 Comprehensive logging**: Structured JSON and text logging with detailed execution tracking and performance monitoring
- **🧹 Smart cleanup**: Configurable cleanup options with debug mode preservation

## 🚀 Quick Start

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/alermol/ChoCallate.git
cd ChoCallate

# Set up the Conda environment
conda env create -f environment.yaml
conda activate ChoCallate
```

### 2. Test Run

```bash
# Run the pipeline on test data
bash run_test.sh

# Optional: Clean up test output
bash cleanup.sh
```

**Note**: The test script expects test data in the `test_data/` directory with the following files:
- `arth_chr1.fasta.gz` — Reference genome (compressed FASTA)
- `test_reads_R1.fq.gz` — Paired-end read 1 (FASTQ)
- `test_reads_R2.fq.gz` — Paired-end read 2 (FASTQ)
- `test_reads_SE.fq.gz` — Single-end reads (FASTQ)
- `sample1.bam` — Example BAM file for BAM input mode

The directory structure should look like:

### 3. Basic Usage

```bash
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv
```

BAM input example (no Bowtie2 index required):

```bash
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --input_format bam \
    --reads_type pe \
    --samples_tsv /path/to/samples_bam.tsv
```

## 🏗️ Pipeline Architecture

### Supported Variant Callers

| Caller      | Diploid Support | Polyploid Support |
|-------------|----------------|-------------------|
| **bcftools**| ✅             | ❌                |
| **GATK4**   | ✅             | ✅                |
| **FreeBayes**| ✅            | ✅                |
| **SNVer**   | ✅             | ✅                |
| **VarDict** | ✅             | ❌                |

### Workflow Scheme

![ChoCallate Pipeline Scheme](ChoCallate_scheme.png)


1. **🔍 Alignment**: Bowtie2-based read alignment with quality filtering and BAM preparation
2. **📊 Coverage Analysis**: Generate coverage information for targeted variant calling
3. **🎯 Zero VCF Generation**: Create baseline VCF with all covered positions
4. **📊 Variant Calling**: Parallel execution of selected variant callers
5. **🤝 Consensus Generation**: Merges results using configurable consensus rules with Python-based SQLite processing
6. **📤 Output**: Final compressed VCF files for SNPs and INDELs

## ⚙️ Configuration

### Essential Parameters

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--reference_genome` | ✅ | - | Reference genome in FASTA format (supports gzipped) |
| `--reference_index` | ✅ for `input_format=fastq` | - | Bowtie2 index prefix for the reference genome (not required for BAM input) |
| `--samples_tsv` | ✅ | `input.tsv` | TSV file with sample information |
| `--input_format` | ✅ | `fastq` | Input files format: `fastq` or `bam` |

### Input/Output Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `ChoCallate_output` | Output directory for results |

### Quality and Filtering Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_coverage` | `5` | Minimum position coverage depth for variant calling |
| `--min_base_quality` | `20` | Minimum base quality for variant calling |
| `--min_map_qual` | `10` | Minimum mapping quality for read filtering |
| `--min_snp_qual` | `20` | Minimum variant quality threshold |

### Data Type Parameters

| Parameter | Default | Choices | Description |
|-----------|---------|---------|-------------|
| `--input_format` | `fastq` | `fastq`, `bam` | Selects whether samples TSV lists FASTQ reads or BAM files |
| `--reads_type` | `pe` | `pe`, `se`, `mx` | Read type: paired-end, single-end, or mixed |
| `--reads_source` | `gbs` | `gbs`, `wgs` | Data source: GBS or whole genome sequencing |
| `--ploidy` | `2` | `≥2` | Ploidy level of the organism |

### Variant Calling Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--effective_callers` | `-` | Comma-separated list of variant callers to use. Use `-` for automatic selection based on ploidy |
| `--cons_type` | `mj` | Consensus type: `mj` (majority), `n1` (n-1), `fc` (full consensus) |

### Resource Allocation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bowtie2_cpu` | `10` | Number of threads for Bowtie2 alignment |
| `--bowtie2_forks` | `1` | Number of parallel Bowtie2 processes |
| `--calling_forks` | `1` | Number of parallel variant calling processes |
| `--zero_vcf_cpu` | `1` | Number of threads for zero VCF generation |
| `--zero_vcf_forks` | `1` | Number of parallel zero VCF processes |
| `--cons_cpus` | `5` | Number of threads for consensus generation |
| `--cons_forks` | `1` | Number of parallel consensus processes |
| `--bcftools_cpu` | `1` | Number of threads for bcftools |
| `--vardict_cpu` | `1` | Number of threads for VarDict |

### Processing Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--win_size` | `1000000` | Window size (in bp) for parallel consensus generation |
| `--debug` | `false` | Keep working directory after pipeline completion |

### Cleanup Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--enable_sample_cleanup` | `true` | Enable/disable sample-specific cleanup (false in debug mode) |
| `--cleanup_intermediate_bam` | `true` | Remove intermediate BAM files (false in debug mode) |
| `--cleanup_intermediate_vcf` | `true` | Remove intermediate VCF files (false in debug mode) |
| `--cleanup_intermediate_subfolders` | `true` | Remove intermediate subfolders (false in debug mode) |
| `--cleanup_input_symlinks` | `true` | Remove symlinks to input files (false in debug mode) |

**Note**: The actual default values are dynamically set based on debug mode. When `--debug` is false (production mode), cleanup is enabled. When `--debug` is true, cleanup is disabled to preserve intermediate files for analysis.

### Logging Parameters

| Parameter | Default | Choices | Description |
|-----------|---------|---------|-------------|
| `--log_level` | `INFO` | `DEBUG`, `INFO`, `WARN`, `ERROR`, `FATAL` | Logging level for pipeline execution |
| `--log_format` | `json` | `json`, `text`, `both` | Log output format |
| `--log_timestamp` | `true` | `true`, `false` | Include timestamps in logs |
| `--log_process` | `true` | `true`, `false` | Include process names in logs |
| `--log_sample` | `true` | `true`, `false` | Include sample IDs in logs |
| `--log_file` | `ChoCallate.log` | - | Main log file path |
| `--log_error_file` | `ChoCallate_errors.log` | - | Error log file path |

### Consensus Types

- **`mj` (Majority Rule)**: Variant is called if majority of callers identify it
- **`n1` (N-1 Consensus)**: Variant is called if n-1 callers identify it (where n is total number of callers)
- **`fc` (Full Consensus)**: Variant is called only if all callers identify it

### Consensus Implementation

The consensus generation uses a sophisticated approach:
- **Zero VCF Integration**: All covered positions from the zero VCF are included in the final output
- **SQLite Processing**: Python scripts use SQLite databases for efficient variant comparison and consensus calculation
- **Window-based Processing**: Genomic regions are processed in parallel using configurable window sizes
- **Quality Filtering**: Variants are filtered based on quality scores and caller agreement

### Automatic Caller Selection

When `--effective_callers` is set to `-` (default), ChoCallate automatically selects appropriate callers among available:
- **Diploid (ploidy=2)**: Uses `bcftools,gatk,freebayes,snver,vardict`
- **Polyploid (ploidy>2)**: Uses `gatk,freebayes,snver` (polyploid-compatible callers only)

### Ploidy and Caller Selection Examples

```bash
# Diploid species (default)
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --ploidy 2 \
    --cons_type mj

# Polyploid species
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --ploidy 4 \
    --cons_type n1 \
    --effective_callers gatk,freebayes,snver
```

## 📁 Input Data Structure

### Samples TSV Format

The structure of `--samples_tsv` depends on two parameters:
- `--input_format`: `fastq` (raw reads) or `bam` (pre-aligned)
- `--reads_type`: `pe`, `se`, or `mx` (applies to FASTQ mode)

Notes (applies to all modes):
- No header line is expected; do not include a header row.
- Fields must be separated by a single TAB character (TSV), not spaces or commas.

#### FASTQ mode (`--input_format fastq`)

Provide 4 columns per sample: `sample_id`, `R1`, `R2`, `SE`.
- Required columns by `--reads_type`:
  - `pe`: columns 1,2,3 required (R1,R2). Column 4 can be `-`.
  - `se`: columns 1 and 4 required (SE). Columns 2 and 3 can be `-`.
  - `mx`: all 4 columns required (R1,R2,SE).

Examples:
```bash
# reads_type=pe
sample1    /path/R1.fq.gz    /path/R2.fq.gz    -

# reads_type=se
sample2    -                 -                 /path/SE.fq.gz

# reads_type=mx
sample3    /path/R1.fq.gz    /path/R2.fq.gz    /path/SE.fq.gz
```

Accepted read formats: `.fq.gz`, `.fastq.gz`, `.fq`, `.fastq`.

#### BAM mode (`--input_format bam`)

Provide at least 2 columns per sample: `sample_id`, `bam_path`. Columns 3–4 are ignored.
- `--reads_type` is accepted but does not affect BAM mode parsing.

Example:
```bash
sample1    /abs/path/sample1.bam    x    x
```

Notes:
- Column 2 must be a valid `.bam` file.

#### File Format Support
- Input reads (FASTQ mode): `.fq.gz`, `.fastq.gz`, `.fq`, `.fastq`
- Reference genome: `.fasta`, `.fa`, `.fna` (compressed or uncompressed)
- Output VCFs: `.vcf.gz` (bgzipped)

### Reference Requirements

- **Format**: FASTA (supports both compressed and uncompressed)
- **Index**: Pre-built Bowtie2 index (required only for `--input_format fastq`)
- **Path**: Absolute paths required

## 📊 Output Structure

```
ChoCallate_output/
├── sample1/
│   ├── sample1.snps.vcf.gz      # Final SNPs VCF (bgzipped)
│   └── sample1.indels.vcf.gz    # Final INDELs VCF (bgzipped)
├── sample2/
│   ├── sample2.snps.vcf.gz
│   └── sample2.indels.vcf.gz
├── ChoCallate_errors.log         # Error log for the entire pipeline
├── ChoCallate.log                # Main log file for the pipeline
├── pipeline_report.html          # Pipeline summary report (HTML)
├── timeline_report.html          # Timeline of process execution (HTML)
└── trace.txt                     # Detailed process trace file
```

## 🔧 Advanced Configuration

### Quality Filtering

```bash
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --min_coverage 10 \
    --min_base_quality 30 \
    --min_map_qual 20 \
    --min_snp_qual 30
```

### Resource Allocation

```bash
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --bowtie2_cpu 16 \
    --cons_cpus 8 \
    --win_size 2000000
```

### Memory Optimization

For large genomes or high read counts, adjust memory allocation:

```bash
# Replace N with desired RAM in GB
sed -i 's/-Xmx1g/-XmxNg/' $CONDA_PREFIX/bin/snver
sed -i 's/-Xmx8g/-XmxNg/' $CONDA_PREFIX/bin/vardict-java
```

## 🏗️ Project Structure

```
ChoCallate/
├── main.nf                      # Main Nextflow pipeline script
├── nextflow.config              # Pipeline configuration
├── functions/                   # Utility functions
│   ├── utils.nf                 # Parameter validation functions
│   └── logging.nf               # Logging utilities
├── bin/                         # Pipeline scripts and variant caller wrappers
│   ├── bcftools_caller.sh       # BCFtools variant calling
│   ├── gatk4_caller.sh          # GATK4 variant calling
│   ├── freebayes_caller.sh      # FreeBayes variant calling
│   ├── snver_caller.sh          # SNVer variant calling
│   ├── vardict_caller.sh        # VarDict variant calling
│   ├── consensus_generation.sh  # Consensus generation script
│   ├── prepare_bam.sh           # BAM preparation and alignment script
│   ├── process_snps.py          # Python script for SNPs consensus
│   ├── process_indels.py        # Python script for indels consensus
│   └── logging_utils.py         # Python logging utilities
├── run_test.sh                  # Test execution script
├── cleanup.sh                   # Test cleanup script
└── README.md                    # This file
```

## 🛠️ Dependencies

All dependencies are managed via Conda:

```yaml
# Core variant callers
- freebayes>=1.3.9
- gatk4=4.6.*
- snver=0.5.3
- vardict-java=1.8.3
- bcftools>=1.20

# Alignment and processing
- bowtie2
- samtools>=1.21
- bedtools
- bedops>=2.4.42

# Pipeline framework
- nextflow
- python
- tabix>=1.11
- parallel
```

## 📈 Performance and Optimization

### Parallel Processing

- **Parallel processing**: Efficient parallel execution of variant callers
- **Multi-threading**: Configurable CPU allocation per process
- **Fork-based parallelism**: Multiple parallel instances for I/O intensive tasks
- **Window-based consensus**: Parallel processing of genomic regions

### Memory Management

- **Streaming consensus**: Efficient memory usage for large datasets
- **Intermediate cleanup**: Automatic removal of temporary files
- **Memory optimization**: Efficient memory usage for large datasets
- **SQLite-based consensus**: Fast in-memory database operations using Python scripts
- **Window-based processing**: Parallel processing of genomic regions for memory efficiency

## 🐛 Troubleshooting

### Common Issues

1. **Memory errors**: Increase memory allocation for SNVer/VarDict
2. **Disk space**: Monitor available disk space for intermediate files
3. **Reference format**: FASTA can be compressed or uncompressed
4. **Path issues**: Use absolute paths for input files

### Debug Mode

```bash
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --debug \
    --log_level DEBUG
```

Debug mode preserves all intermediate files for analysis.

### Cleanup Options

```bash
# Disable cleanup for debugging
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --enable_sample_cleanup false \
    --debug

# Custom cleanup configuration
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv samples.tsv \
    --cleanup_intermediate_bam false \
    --cleanup_intermediate_vcf true
```


## 📄 Citation

**APA Style**:  
Ermolaev, A. (2025). *ChoCallate: Consensus variant calling pipeline* [Computer software]. GitHub. https://github.com/alermol/ChoCallate

**BibTeX**:  
```bibtex
@software{ChoCallate,
  author = {Ermolaev, A.},
  title = {ChoCallate: Consensus variant calling pipeline},
  url = {https://github.com/alermol/ChoCallate},
  year = {2025}
}
```

## 🗺️ Development Roadmap

ChoCallate is actively developed with a clear vision for future enhancements. Here's our roadmap for upcoming versions:

- Add New Germline Variant Callers
- Add New Short Read Mapping Tools
- Add Somatic Variant Callers
- Add Long-Read Variant Callers
- Add Long-Read Mapping Tools
- Add AI-Powered Features
    - ML-based automatic consensus generation
    - AI-powered variant quality assessment
- Add Containerized Solution


### 🛠️ Development Priorities

1. **Performance Optimization**: Implement advanced strategies to significantly reduce pipeline runtime
2. **Error Handling**: Improved error recovery and user feedback
3. **New Variant Callers**: Integration of cutting-edge tools
4. **Quality Metrics**: Enhanced quality assessment and reporting
5. **Format Support**: Additional input/output format compatibility
6. **Testing**: Comprehensive test suite and validation

### 🤝 Contributing to Development

We welcome contributions from the community! Here's how you can help:

#### Development Areas
- **Core Pipeline**: Nextflow workflow optimization
- **Variant Callers**: Integration of new variant calling tools
- **Consensus Algorithms**: Improved consensus generation methods
- **Quality Control**: Enhanced quality assessment tools
- **Documentation**: User guides and technical documentation

#### Getting Started
1. Fork the repository
2. Create a feature branch
3. Implement your changes
4. Add tests and documentation
5. Submit a pull request

---

## 📜 License

MIT License - see [LICENSE](LICENSE) file for details.

---

**Need help?** Open an issue on GitHub or check our troubleshooting guide above.
