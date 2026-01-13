# ChoCallate 🍫

**ChoCallate** (**Cho**rus of **Call**ers) is a high-performance, automated pipeline for consensus-based variant calling that combines multiple variant callers to produce robust, high-confidence single-nucleotide variants (SNVs) and indels (INDELs).

## What is ChoCallate?

ChoCallate addresses a critical challenge in variant calling: individual variant callers can produce different results for the same genomic data, leading to uncertainty in variant identification. By implementing a **consensus-driven approach**, ChoCallate combines results from multiple state-of-the-art variant callers and applies configurable consensus rules to generate reliable, high-quality variant calls.

### Key Features

- **Consensus-driven approach**: Combines multiple variant callers using configurable consensus rules
- **Ploidy flexibility**: Supports both diploid and polyploid species with automatic caller selection
- **Multiple consensus types**: Majority rule, n-1 consensus, and full consensus options
- **Dual input support**: Processes both FASTQ (raw reads) and BAM (pre-aligned) files, allowing flexible integration of sequencing data at different analysis stages
- **Flexible input compatibility**: Works with GBS (Genotyping-by-Sequencing) and WGS data
- **Parallel processing**: Efficient parallel execution for optimal performance
- **Configurable quality filtering**: Multiple filtering steps based on coverage, base quality, and SNP quality
- **Comprehensive logging**: Structured JSON and text logging with detailed execution tracking and performance monitoring
- **Smart cleanup**: Configurable cleanup options with debug mode preservation
- **BCF-native processing**: Uses compressed BCF format throughout the pipeline for optimal performance (optional VCF output available)
- **Optional single-file output**: Merge per-sample results into single multi-sample BCF

## Quick Start

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
    --samples_tsv /path/to/samples_bam.tsv
```

### 4. Command-Line Help

Get help and version information:

```bash
# Show version information
nextflow run main.nf --version

# Show help
nextflow run main.nf --help
```

## Pipeline Architecture

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


1. **Alignment**: Bowtie2-based read alignment with quality filtering and BAM preparation
2. **Coverage Analysis**: Generate coverage information for targeted variant calling
3. **Zero BCF Generation**: Create position-template (zero) BCF with all covered positions
4. **Variant Calling**: Parallel execution of selected variant callers
5. **Consensus Generation**: Merges results using configurable consensus rules with Python-based SQLite processing
6. **Optional Merge Step**: When enabled, merges all samples' BCFs into single final SNP/INDEL BCFs
7. **Output**: Final compressed BCF (or VCF when `--output_vcf` is enabled) files for SNPs and INDELs

## Configuration

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
| `--output_vcf` | `false` | Output compressed VCF instead of compressed BCF |

### Quality and Filtering Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_coverage` | `5` | Minimum position coverage depth for variant calling |
| `--min_base_quality` | `5` | Minimum base quality for variant calling |
| `--min_map_qual` | `5` | Minimum mapping quality for read filtering |
| `--min_snp_qual` | `5` | Minimum variant quality threshold |

### Data Type Parameters

| Parameter | Default | Choices | Description |
|-----------|---------|---------|-------------|
| `--input_format` | `fastq` | `fastq`, `bam` | Selects whether samples TSV lists FASTQ reads or BAM files |
| `--reads_source` | `gbs` | `gbs`, `wgs` | Data source: GBS or whole genome sequencing |
| `--ploidy` | `2` | `≥2` | Ploidy level of the organism |

### Variant Calling Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--effective_callers` | `-` | Comma-separated list of variant callers to use (case-insensitive). Use `-` for automatic selection based on ploidy. |
| `--cons_type` | `mj` | Consensus type: `mj` (majority), `n1` (n-1), `fc` (full consensus) |

### Resource Allocation Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--bowtie2_cpu` | `10` | Number of threads for Bowtie2 alignment |
| `--bowtie2_forks` | `1` | Number of parallel Bowtie2 processes |
| `--calling_forks` | `1` | Number of parallel variant calling processes |
| `--zero_bcf_cpu` | `1` | Number of threads for zero BCF generation |
| `--zero_bcf_forks` | `1` | Number of parallel zero BCF processes |
| `--cons_cpus` | `5` | Number of threads for consensus generation |
| `--cons_forks` | `1` | Number of parallel consensus processes |
| `--bcftools_cpu` | `1` | Number of threads for bcftools |
| `--vardict_cpu` | `1` | Number of threads for VarDict |
| `--merge_bcfs_cpus` | `1` | Number of threads for BCF merge step |

### Processing Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--win_size` | `1000000` | Window size (in bp) for parallel consensus generation |
| `--test_run` | `false` | Enable test mode: process only the first `test_run_limit` samples from `--samples_tsv` |
| `--test_run_limit` | `2` | Maximum number of samples to process when `--test_run` is enabled |
| `--debug` | `false` | Keep working directory after pipeline completion |
| `--bowtie2_extra_args` | `""` | Extra arguments passed directly to Bowtie2 during alignment (use as is) |
| `--gatk_leftalignindels_extra_args` | `""` | Extra arguments passed to gatk LeftAlignIndels (use as is) |
| `--bcftools_mpileup_extra_args` | `""` | Extra arguments appended to `bcftools mpileup` (use as is)|
| `--bcftools_call_extra_args` | `""` | Extra arguments appended to `bcftools call` (use as is)|
| `--freebayes_extra_args` | `""` | Extra arguments appended to `freebayes` (use as is)|
| `--gatk4_extra_args` | `""` | Extra arguments appended to `gatk HaplotypeCaller` (use as is)|
| `--snver_extra_args` | `""` | Extra arguments appended to `snver` (use as is)|
| `--vardict_extra_args` | `""` | Extra arguments appended to `vardict-java` (use as is)|
| `--bcftools_merge_extra_args` | `""` | Extra arguments appended to `bcftools merge` (use as is)|
| `--merge_bcfs_forks` | `1` | Number of parallel merge processes |
| `--single_file` | `false` | If `true`, output one merged pair of final BCFs |

### Cleanup Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--enable_sample_cleanup` | `true` | Enable/disable sample-specific cleanup (false in debug mode) |
| `--cleanup_intermediate_bam` | `true` | Remove intermediate BAM files (false in debug mode) |
| `--cleanup_intermediate_bcf` | `true` | Remove intermediate BCF files (false in debug mode) |
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

### Help and Version Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--help` | `false` | Show help message and exit |
| `--version` | `false` | Show version information and exit |

### Consensus Types

- **`mj` (Majority Rule)**: Variant is called if majority of callers identify it
- **`n1` (N-1 Consensus)**: Variant is called if n-1 callers identify it (where n is total number of callers)
- **`fc` (Full Consensus)**: Variant is called only if all callers identify it

### Consensus Implementation

The consensus generation uses a sophisticated approach:
- **Zero BCF Integration**: All covered positions from the zero BCF are included in the final output
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

## Input Data Structure

### Samples TSV Format

The structure of `--samples_tsv` depends on the `--input_format` parameter:
- `--input_format`: `fastq` (raw reads) or `bam` (pre-aligned)

Notes (applies to all modes):
- No header line is expected; do not include a header row.
- Fields must be separated by a single TAB character (TSV), not spaces or commas.

#### FASTQ mode (`--input_format fastq`)

Provide 4 columns per sample: `sample_id`, `R1`, `R2`, `SE`.
- At least one of columns 2, 3, or 4 must contain valid FASTQ file paths
- **Valid read combinations:**
  - **R1 + R2 only**: Paired-end reads (column 4 can be `-`)
  - **SE only**: Single-end reads (columns 2 and 3 can be `-`)
  - **R1 + R2 + SE**: Mixed reads (all three types)
- If both R1 and R2 are provided, they must have equal counts
- Empty columns can be marked with `-`

Examples:
```bash
# Paired-end reads
sample1    /path/R1.fq.gz    /path/R2.fq.gz    -

# Single-end reads
sample2    -                 -                 /path/SE.fq.gz

# Mixed reads (paired-end + single-end)
sample3    /path/R1.fq.gz    /path/R2.fq.gz    /path/SE.fq.gz
```

Accepted read formats: `.fq.gz`, `.fastq.gz`, `.fq`, `.fastq`.

#### BAM mode (`--input_format bam`)

Provide at least 2 columns per sample: `sample_id`, `bam_path`. Columns 3–4 are ignored.

Example:
```bash
sample1    /abs/path/sample1.bam    -    -
```

Notes:
- Column 2 must be a valid `.bam` file.

#### File Format Support
- Input reads (FASTQ mode): `.fq.gz`, `.fastq.gz`, `.fq`, `.fastq`
- Reference genome: `.fasta`, `.fa`, `.fna` (gzipped or ungzipped)
- Variant caller output: `.bcf` (compressed BCF format)
- Final output: `.bcf` (compressed BCF format) or `.vcf.gz` (compressed VCF format when `--output_vcf` is enabled)

### Reference Requirements

- **Format**: FASTA (supports both compressed and uncompressed)
- **Index**: Pre-built Bowtie2 index (required only for `--input_format fastq`)
- **Path**: Absolute paths required

## Output Structure

Default (per-sample outputs):

```
ChoCallate_output/
├── sample1/
│   ├── sample1.snps.bcf      # Final SNPs BCF (compressed)
│   └── sample1.indels.bcf    # Final INDELs BCF (compressed)
├── sample2/
│   ├── sample2.snps.bcf
│   └── sample2.indels.bcf
├── ChoCallate_errors.log         # Error log for the entire pipeline
├── ChoCallate.log                # Main log file for the pipeline
├── pipeline_report.html          # Pipeline summary report (HTML)
├── timeline_report.html          # Timeline of process execution (HTML)
└── trace.txt                     # Detailed process trace file
```

Single-file mode (`--single_file`):

```
ChoCallate_output/
├── final.snps.bcf               # Merged SNPs across all samples
├── final.indels.bcf             # Merged INDELs across all samples
├── ChoCallate_errors.log
├── ChoCallate.log
├── pipeline_report.html
├── timeline_report.html
└── trace.txt
```

## Advanced Configuration

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

### Test Run Mode

Use test mode to quickly validate configuration on a small subset of samples:

```bash
nextflow run main.nf \
    --reference_genome /path/to/reference.fasta \
    --reference_index /path/to/reference_index \
    --samples_tsv /path/to/samples.tsv \
    --test_run true \
    --test_run_limit 2
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

## Project Structure

```
ChoCallate/
├── main.nf                      # Main Nextflow pipeline script
├── nextflow.config              # Pipeline configuration
├── environment.yaml             # Conda environment specification
├── LICENSE                      # MIT License file
├── functions/                   # Utility functions
│   ├── utils.nf                 # Parameter validation functions
│   ├── logging.nf               # Logging utilities
│   ├── help_version.nf          # Help and version display module
│   ├── calling.nf               # Variant calling workflow
│   ├── prepare_bam.nf           # BAM preparation workflow
│   ├── coverage_generation.nf   # Coverage analysis workflow
│   ├── create_fai_index.nf      # FASTA index creation
│   ├── create_seq_dict.nf       # Sequence dictionary creation
│   ├── generate_zero_bcf.nf     # Zero BCF generation workflow
│   ├── generate_consensus.nf    # Consensus generation workflow
│   ├── merge_bcfs.nf            # Merge per-sample BCFs into single outputs
│   └── cleanup_sample_temp.nf   # Sample cleanup workflow
├── bin/                         # Pipeline scripts and variant caller wrappers
│   ├── bcftools_caller.sh       # BCFtools variant calling
│   ├── gatk4_caller.sh          # GATK4 variant calling
│   ├── freebayes_caller.sh      # FreeBayes variant calling
│   ├── snver_caller.sh          # SNVer variant calling
│   ├── vardict_caller.sh        # VarDict variant calling
│   ├── consensus_generation.sh  # Consensus generation script
│   ├── prepare_bam.sh           # BAM preparation and alignment script
│   ├── process_snps.py          # Python script for SNPs consensus
│   └── process_indels.py        # Python script for indels consensus
├── run_test.sh                  # Test execution script
├── cleanup.sh                   # Test cleanup script
└── README.md                    # This file
```

## Dependencies

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

## Troubleshooting

### Common Issues

1. **Memory errors**: Increase memory allocation for SNVer/VarDict
2. **Disk space**: Monitor available disk space for intermediate files
3. **Path issues**: Use absolute paths for input files

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
    --cleanup_intermediate_bcf true
```


## Citation

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

## Development Roadmap

ChoCallate is actively developed with a clear vision for future enhancements. Here's a roadmap for upcoming versions:

- Add New Germline Variant Callers
- Add New Short Read Mapping Tools
- Add Somatic Variant Callers
- Add Long-Read Variant Callers
- Add Long-Read Mapping Tools
- Add AI-Powered Features
    - ML-based automatic consensus generation
    - AI-powered variant quality assessment
- Add Containerized Solution


### Development Priorities

1. **Performance Optimization**: Implement advanced strategies to significantly reduce pipeline runtime
2. **Error Handling**: Improved error recovery and user feedback
3. **New Variant Callers**: Integration of cutting-edge tools
4. **Quality Metrics**: Enhanced quality assessment and reporting
5. **Format Support**: Additional input/output format compatibility

### Contributing to Development

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
4. Add documentation
5. Submit a pull request

---

## License

MIT License - see [LICENSE](LICENSE) file for details.

---

**Need help?** Open an issue on GitHub or check our troubleshooting guide above.
