# ChoCallate 🍫

[![GitHub Release](https://img.shields.io/github/v/release/alermol/ChoCallate)](https://github.com/alermol/ChoCallate/releases) [![Changelog](https://img.shields.io/badge/Changelog-orange)](CHANGELOG.md) [![GitHub License](https://img.shields.io/github/license/alermol/chocallate)](LICENSE) [![Wiki](https://img.shields.io/badge/Wiki-red)](https://github.com/alermol/ChoCallate/wiki)

**ChoCallate** (**Cho**rus of **Call**ers) is a Nextflow pipeline for **consensus-based variant calling**. It maps FASTQ reads or accepts aligned BAMs, calls variants with bcftools, FreeBayes and GATK, and applies a configurable consensus rule to produce high-confidence SNVs and indels.

## Quick start

Requirements: Linux, Conda/Mamba and Git. Nextflow and all analysis tools are installed by the environment. The first run resolves the `nf-schema` and `nf-boost` plugins; an offline Conda setup must install them once, while the Docker image bundles them.

```bash
git clone --depth 1 https://github.com/alermol/ChoCallate.git
cd ChoCallate
conda env create -y -f environment.yaml --no-default-packages
conda activate ChoCallate
cp assets/templates/config.yaml my_run.yaml
```

Minimal `my_run.yaml`:

```yaml
input:
  samples_tsv: /absolute/path/samples.tsv
  reference_genome: /absolute/path/reference.fa.gz
  reference_index_dir: /absolute/path/indexes

output:
  directory: ChoCallate_output
```

`samples.tsv` is tab-separated:

```text
sample1	/path/sample1_R1.fq.gz	/path/sample1_R2.fq.gz
```

Supported layouts are paired-end (`sample_id`, R1, R2), single-end (`sample_id`, R1), mixed (`sample_id`, R1, R2, unpaired; Bowtie2 only) and BAM (`sample_id`, BAM). `reference_index_dir` is required for FASTQ input and unused for BAM input. Index files must match the selected mapper.

Run:

```bash
nextflow run main.nf -params-file my_run.yaml
```

## Docker

```bash
docker pull alermol/chocallate:latest
docker run --rm \
  -u "$(id -u):$(id -g)" \
  -v "${PWD}/input_data:/workspace" \
  -w /workspace \
  alermol/chocallate:latest \
  -params-file config.yaml
```

The image bundles the Nextflow plugins and works offline and under an arbitrary UID. Paths in `config.yaml` and `samples.tsv` must refer to files inside the mounted `/workspace`.

## Outputs

By default, a run publishes one merged multi-sample `consensus.bcf`. Set `output.type: "sample"` for one file per sample, and `output.format: "vcf"` for bgzip-compressed VCF. All covered sites are kept by default, including monomorphic positions; set `output.remove_invariant: true` to drop them. Mapping statistics, `dropped_samples.txt` and standard Nextflow reports are also published.

## Documentation

Detailed documentation is in the [Wiki](https://github.com/alermol/ChoCallate/wiki):

- [Installation](https://github.com/alermol/ChoCallate/wiki/Installing-ChoCallate), including Docker, Apptainer and offline use.
- [Pipeline overview](https://github.com/alermol/ChoCallate/wiki/Pipeline-Overview) and [inputs and reference](https://github.com/alermol/ChoCallate/wiki/Inputs-and-Reference) requirements.
- [CLI reference](https://github.com/alermol/ChoCallate/wiki/CLI-Reference): every parameter, including resources, callers and consensus options.
- [Consensus model](https://github.com/alermol/ChoCallate/wiki/Consensus-Model): voting, caller normalization and degraded mode.
- [Output structure](https://github.com/alermol/ChoCallate/wiki/Output-Structure): files, header metadata and FORMAT/INFO annotations.
- [Merge strategies and testing](https://github.com/alermol/ChoCallate/wiki/Merge-Strategies-and-Testing): multi-sample merge and open-file limits.
- [Operations and troubleshooting](https://github.com/alermol/ChoCallate/wiki/Operations-and-Troubleshooting): resources, retries and common errors.
- [Testing and CI](https://github.com/alermol/ChoCallate/wiki/Testing): pytest, nf-test, snapshots and lint commands.
- [Developer guide](https://github.com/alermol/ChoCallate/wiki/Developer-Guide): repository layout and process conventions.

The full configuration template is [`assets/templates/config.yaml`](assets/templates/config.yaml).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

[MIT](LICENSE)

## Development roadmap

See [ROADMAP.md](ROADMAP.md).
