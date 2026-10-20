# Changelog

All notable changes to this project will be documented in this file.

## [4.0.0] - 2026-10-20

### Added

- `nextflow_schema.json` and the `nf-schema` plugin validate every parameter before the first task: enum values, numeric minimums, input-file existence and unknown keys abort with the offending key and the accepted values. Every parameter now has a default in `nextflow.config` and documentation in `assets/templates/config.yaml`.
- Samplesheet validation: blank lines and comments are skipped, only the first non-comment line can be a header, sample ids must be unique and shell-safe, every row must have the exact column count for `input.format`/`input.reads_type`, and relative read paths are resolved from the samplesheet directory and must exist.
- `calling.<caller>.min_snp_qual` and `calling.<caller>.min_indel_qual` provide per-caller, per-variant-type quality filters and are required for every configured caller.
- `calling.allow_caller_failure` (default `false`) enables degraded mode: a failed caller is skipped with a warning, `CallersFailed` is recorded, the run continues while `consensus.threshold` is satisfiable, and samples with no surviving caller are listed in `dropped_samples.txt` (and in `SamplesDropped` in a merged output).
- Caller `extra_args` that request non-diploid ploidy, a copy-number map or custom ploidy regions are rejected, because consensus genotypes are diploid.
- `consensus.absent_vote` (`ref`/`abstain`), `consensus.filtered_vote` (`abstain`/`ref`) and `consensus.tie_policy` (`nocall`/`first`) configure the voting model.
- `output.compress_reference_blocks` (default `true`) collapses runs of hom-reference positions into multi-base gVCF blocks with `INFO/END` when reference records are written, in both per-sample and merged output.
- `output.publish_mode` (`copy` or `link`, default `copy`) copies or links results instead of moving them out of the work directory; `copy` keeps results independent of work-directory cleanup, while `link` requires the work directory to be retained.
- `input.reference_index_prefix` sets the mapper index prefix when it differs from the reference file name; it must be a plain file name without path separators.
- `mapping.platform` (default `ILLUMINA`) and `mapping.library` (default the sample id) set the read-group `PL` and `LB` values for mapping and BAM input; values are restricted to letters, digits, `_`, `.`, `:` and `-`.
- `mapping.aligner_reserve_mb` reserves memory for the aligner index while `samtools sort` runs in the same task; `0` estimates the reserve from the staged index files plus 1 GiB.
- `resources.low`, `resources.medium` and `resources.high` define per-label CPU, memory and time limits; `resources.limits` caps those values after retry escalation, and JVM heaps are derived from task memory.
- `resources.max_retries` controls retries for transient failures (exit status 137 or 140), with memory and time multiplied by the attempt number.
- `boost.cleanup` is now a params-file option (`"v1"`, `"v2"` or `false`) instead of a hard-coded `nextflow.config` setting.
- `SCATTER_INTERVALS` computes the coverage BED chunks once per sample and shares them between every caller and consensus.
- `EXPORT_OUTPUT` applies `output.format`, `output.remove_invariant` and `output.split_multiallelic` once at the end and publishes per-sample or merged files according to `output.type`; internal consensus and merge files are always BCF, including for `output.format: "vcf"`.
- Mapping statistics are published under `mapping_stats/`: raw `samtools flagstat` plus final `samtools flagstat`, `samtools stats` and `samtools idxstats`, and `mapping_summary.tsv` with per-sample alignment, duplicate, MAPQ, insert-size and error-rate metrics.
- `bin/postprocess_bam.sh` provides shared mapping/BAM post-processing: mate fixing, `samtools sort`, optional `samtools markdup`, MAPQ filtering, optional GATK `LeftAlignIndels` and one CSI index.
- BAM input read groups are normalised with `samtools addreplacerg` (`ID`/`SM` set to the sample id, `PL`/`LB` from `mapping.platform`/`mapping.library`), so GATK and `bcftools reheader` see consistent read groups.
- BAM inputs and mapping output are checked against the reference: `@SQ` names and lengths must match the `.fai`, and when the BAM header carries `M5` checksums they are compared with the sequence dictionary; mismatches abort before calling.
- The generated `.fai` is checked for contigs longer than 2,147,483,647 bp (2^31-1), because BAM/BCF positions and GATK use 32-bit coordinates.
- `bin/filter_bed.py` strips `track`, `browser` and `#` lines before `sort-bed` and rejects contigs absent from the reference `.fai`; consensus rejects coverage BED intervals with negative or reversed coordinates or ends past the contig.
- `coverage.on_empty` (`fail` by default, or `skip`) controls what happens when no base of a sample reaches `coverage.min_coverage`; skipped samples are listed in `dropped_samples.txt` (and in `SamplesDropped` in a merged output).
- `bin/merge_outputs.sh` is the production merge chain used by `MERGE_OUTPUTS` and the pytest merge tests; `bin/drop_refonly_duplicates.py` removes the redundant `<*>` record created by `bcftools merge --gvcf`.
- `bin/merge_cho_headers.py` builds the merged `##ChoCallate=` line from the union of per-sample `Callers` and `CallersFailed`, adds `SamplesDropped` when samples were dropped, and writes `##ChoCallateSamples=` when caller sets differ or any caller failed.
- The merge step checks the process open-file limit before the first task for `output.type: "single"` and aborts with the required `ulimit -n` value when the sample count and `consensus.cpu` would exceed it.
- `bin/build_consensus_accel.py` builds an optional Cython accelerator for `generate_consensus.py`; the Docker image builds it automatically, and an equivalent pure-pysam implementation remains as a fallback.
- One structured `##ChoCallate=` header line per consensus file records the version, threshold, callers, rule, tie policy, absence policy, filtered-record policy and reference-block flag.
- The `tests/python` pytest suite covers consensus (pure and accelerated), merge, BAM post-processing, validation, resources, retries, resume/cache and cleanup; the `tests/nf` nf-test suite contains module tests and end-to-end pipeline cases on the bundled `test_run` data.
- GitHub Actions CI runs `ruff`, `nextflow lint`, `actionlint`, pytest with and without the Cython accelerator, nf-test, and Docker image smoke tests including offline plugin resolution.

### Changed

- `output.remove_invariant` now defaults to `false`, so a default run keeps monomorphic/reference sites in both per-sample and merged BCF/VCF output.
- Hom-reference positions are written as gVCF-style `ALT="<*>"` records with `INFO/END` instead of `ALT="."`; single-base records carry `AD=(reference depth, other-allele depth)`, while multi-base blocks omit `AD` and carry the block minimum in `FORMAT/DP`. `MERGE_OUTPUTS` uses `bcftools merge --gvcf`, so merged `FORMAT/AD` reports honest ALT depth for hom-reference samples.
- `coverage.min_coverage` now defaults to `5` (was `2`), must be at least `1`, and `bam_filter.min_map_qual` now defaults to `10` (was `0`).
- `rmdup.mode` (`none`, `mark` or `remove`; default `mark`) replaces `rmdup.enabled` and `picard MarkDuplicates`; `samtools markdup` flags duplicates by default, and callers, coverage and consensus skip duplicate-flagged reads.
- Reference compression is detected from the BGZF header, so `ref_genome.bgzip` is no longer needed; plain gzip references are rejected while plain and BGZF-compressed FASTA are accepted.
- `input.reference_index_dir` is required only for `format: "fastq"`, so BAM-only runs no longer need mapper indexes.
- Mapper indexes are resolved from `input.reference_index_dir` and staged as task inputs instead of being read from that host path at task runtime, so mapping works in containers, on cloud executors and with `scratch`; an incomplete or missing index aborts before the first task.
- Consensus inputs are staged as `tmp/callers/<caller>.bcf` and passed as `--input caller=path`, so caller identity and voting order no longer depend on positional file order when a caller fails.
- Callers are normalised with `bcftools norm -f ref -m -any --atomize`; multi-allelic and complex records are split, low-QUAL records are kept with `FILTER=LowQual`, and the old `--min-alleles 2 --max-alleles 2` filter is no longer applied.
- Consensus is now two-stage and caller-order independent by default: an ALT becomes a candidate once `consensus.threshold` callers support it, then the majority zygosity among its supporters decides the genotype; ties are a no-call unless `consensus.tie_policy: "first"`, and multiple supported ALT alleles produce ploidy-consistent split records.
- Filtered caller records abstain from the vote by default (`consensus.filtered_vote: "abstain"`), no-calls no longer count as hom-reference votes, and positions with an IUPAC-ambiguous reference base or caller allele are omitted from consensus output.
- `FORMAT/DP` counts collapsed observations (one per read name) instead of every covering read, matching the overlap-corrected coverage BED; a deletion is not an allele, so an overlapping mate's base is preferred when tallying reference `AD`.
- `calling.bcftools.mpileup.max_depth` and `max_idepth` (default `100000`) override the bcftools default of 250; FreeBayes receives `--use-best-n-alleles` (default `4`) and `--limit-coverage` (default `100000`); GATK receives `--standard-min-confidence-threshold-for-calling min(min_snp_qual, min_indel_qual)`.
- GATK calling intervals are padded by `calling.gatk.interval_padding` (default `150`), and each padded chunk's calls are restricted to its own unpadded BED interval before concatenation.
- Coverage is computed with `mosdepth --no-per-base --quantize` and sorted with `sort-bed`/`bedops` instead of a per-base `samtools depth | awk` stream.
- Mapping and BAM preparation share `bin/postprocess_bam.sh`; in mapping tasks `mapping.cpu` is split between mapper and post-processing, sort memory is derived from task memory after the aligner reserve, and each BAM is indexed once.
- Every process now uses `process_low`, `process_medium` or `process_high`; every JVM receives an `-Xmx` derived from task memory, `Runtime.runtime.availableProcessors()` is no longer used, `maxForks` directives are gone, and caller failures are no longer ignored unconditionally.
- All process scripts run under `bash -euo pipefail`, process-local `$TMPDIR` is used for temporary files, and cleanup removes `tmp/`, `index/` and `TMPDIR` only after a successful task.
- The bundled `test_run` dataset is reduced to a 200-kb slice and 20,000 read pairs, so the end-to-end test takes about 3 minutes.
- `environment.yaml` adds `mosdepth`, `cython`, `setuptools` and `c-compiler` and updates the tool stack; `htslib` provides `bgzip`/`tabix`.
- The Docker image installs the pipeline under `/opt/chocallate`, pins the `nf-schema` and `nf-boost` plugins at build time for offline runs, and provides a world-writable `HOME`, so it works with Docker, Apptainer and Singularity under an arbitrary UID.
- `nextflow.config` now declares defaults for every parameter plus `manifest.name`, `nextflowVersion`, `homePage` and `mainScript`.
- `generate_consensus.py` streams one cursor per caller BCF instead of materialising all records of a coverage span and is about 20x faster on a fragmented 2.1 Mbp coverage BED through span grouping, a fast path for positions without caller records and bulk pileup accessors; an optional Cython accelerator speeds up the pileup and record writing further, with an equivalent pure-pysam fallback.
- `merge_bcf_tree.py` handles a single input file and aborts pending FIFO readers/writers when the final merge fails.

### Removed

- `ref_genome.bgzip`; compression is detected automatically.
- `rmdup.enabled`; use `rmdup.mode` instead.
- The shared `calling.min_snp_qual` key; configure `calling.<caller>.min_snp_qual` and `calling.<caller>.min_indel_qual` instead.
- `bin/fill_missing_ad.py`; the gVCF-style `<*>` representation makes post-merge AD filling unnecessary.
- The standalone `tabix` and unused `bedtools` conda dependencies.
- The `--remove_invariant` and `--split_multiallelic` options of `generate_consensus.py`; `EXPORT_OUTPUT` now applies both transformations.

### Fixed

- `generate_consensus.py` no longer drops the last base of every coverage BED interval, so 1-bp intervals produce output and variants at interval ends are retained.
- Deletion handling no longer skips the deletion REF span or loses overlapping caller records: the walk advances one base at a time, caller records inside a deletion are still voted on, and reference records are not written over deleted bases.
- ALT support is keyed by `(REF, ALT)` instead of ALT alone, so deletions of different lengths and atomised records with different REF alleles at the same position remain distinct and caller-order independent.
- `generate_consensus.py` no longer crashes with `TypeError` on a caller record with a missing genotype (`./.`, `0/.`).
- `MERGE_OUTPUTS` keeps monomorphic sites when `output.remove_invariant: false`, matching the documented meaning; the merged output no longer drops every reference-only record.
- `output.type: "single"` with `output.format: "vcf"` no longer fails by passing per-sample `.vcf.gz` files to a merge stage that expects BCF; internal files are always BCF and the VCF is created only by the export stage.
- `PREPARE_BAM` declares the BAM it actually writes (`mapping.bam`), so `input.format: "bam"` no longer fails with "Missing output file(s)".
- `merge_bcf_tree.py` no longer fails with "Expected two or more files to merge" when it runs with a single chunk (`--cpus 1`), and a failed final merge no longer hangs on FIFOs.
- `bin/stage_cleanup.sh` preserves staged inputs and temporary files when a task fails and never deletes a `TMPDIR` outside the task work directory.
- `GENERATE_COVERAGE` stages `include_bed` and `exclude_bed` in separate directories, so BED inputs with the same file name no longer overwrite each other.
- `test_run/run_test.sh` sets `NXF_OFFLINE` (the variable was misspelled `NFX_OFFLINE`), so the offline setting takes effect.


## [3.2.3] - 2026-10-01

### Fixed

- Each sample's BAM is now joined with its own coverage BED by `sample_id` before variant calling and consensus generation. `GENERATE_COVERAGE` no longer drops `sample_id` from its output, so the coverage BED can no longer be paired with another sample's BAM and callers no longer report variants at positions below `coverage.min_coverage`.
- Merged multi-sample outputs now have a deterministic sample order that follows `samples_tsv`. `MERGE_OUTPUTS` ranks the per-sample consensus files by their `samples_tsv` position instead of relying on process completion order, and the merge list is sorted naturally (`ls -1v`) so the order is preserved for ≥ 10 samples.
- BAM input (`format: "bam"`) now works: the sample channel no longer reads a separate `input.input_format` key that disagreed with the `input.format` key used everywhere else.
- `CLIParamsValidation.effective_callers_validation` now actually validates caller names and no longer calls `split()` before its null check; `cons_threshold_validation` also checks `callers` for null before using it.

### Changed

- The documented input format key is now `format` (was `input_format`) in the README, matching `assets/templates/config.yaml`.

## [3.2.2] - 2026-09-30

### Fixed

- `gatk SplitIntervals` in all calling processes (`CALLING_BCFTOOLS`, `CALLING_FREEBAYES`, `CALLING_GATK`) now writes its temp files to the process `$TMPDIR` (`--java-options "-Djava.io.tmpdir=$TMPDIR"` and `--tmp-dir $TMPDIR`) instead of the system temporary directory.

## [3.2.1] - 2026-09-30

### Changed

- The `bcftools` caller now emits `FORMAT/AD` and `FORMAT/DP`, which are requested directly from `bcftools mpileup` (`-a FORMAT/AD,FORMAT/DP`) instead of being recovered downstream.
- The reference sequence dictionary (`ref_genome.dict`) is now passed to all calling processes and to `GENERATE_CONSENSUS`, as required by `gatk SplitIntervals`.
- All calling processes (`CALLING_BCFTOOLS`, `CALLING_FREEBAYES`, `CALLING_GATK`) and `GENERATE_CONSENSUS` now split the coverage BED with `gatk SplitIntervals` (balanced, size-aware intervals) instead of `split`/`bedops --chop` + `split`.
- `FORMAT/AD` and `FORMAT/DP` recovery in `generate_consensus.py` is faster: the genotype walk resumes from the last record instead of rescanning from the start at every position, and the BAM pileup is streamed position-by-position instead of materialising the whole region (bounded memory, no second pass).

### Fixed

- Removal of monomorphic (invariant) positions now filters hom-ref sites with `COUNT(GT="RR")=N_SAMPLES` instead of counting all-homozygous/all-heterozygous genotypes, so only truly invariant positions are dropped and genuine variants are no longer removed.

## [3.2.0] - 2026-09-04

### Added

- Consensus outputs now include `FORMAT/AD` and `FORMAT/DP`, recovered from the sample BAM pileup (overlapping PE mates are deduplicated for allele counts; `DP` counts all covering reads to match `samtools depth` used for the coverage BED).
- Site-level `INFO` annotations via `bcftools +fill-tags` on both per-sample and merged outputs: `AN`, `AC`, `AF`, `NS`, `AC_Hom`, `AC_Het`, `MAF`, `TYPE`, `F_MISSING`, and `DP` (sum of `FORMAT/DP`).

### Changed

- Consensus `QUAL` is always missing (`.`) instead of a placeholder score.
- `GENERATE_CONSENSUS` takes the sample BAM as input so AD/DP can be recovered during consensus writing.

## [3.1.0] - 2026-09-01

### Added

- Wiki: [Merge Strategies and Testing](https://github.com/alermol/ChoCallate/wiki/Merge-Strategies-and-Testing).

### Changed

- `MERGE_OUTPUTS` now uses chunked pipe merge for N samples ≥ 2 (no intermediate files, named pipes only).

## [3.0.0] - 2026-08-27

### Changed

- Dependency version pins in `environment.yaml` were relaxed to minimum versions (`>=`), and the Nextflow requirement was raised to `>=26.04.6`.
- The merging of single-sample VCF/BCF files into multi-sample files was significantly accelerated by using hierarchical merging.
- `MERGE_OUTPUTS` now uses `consensus.cpu` for thread allocation. The previous retry strategy that decreased CPU count on each attempt was removed.
- Default `output.type` changed from `"sample"` to `"single"` (one multi-sample consensus file).

### Added

- Optional removal of monomorphic (invariant) positions in both single- and multi-sample outputs to reduce file size (`output.remove_invariant`).
- Optional splitting of multi-allelic variants in output multi-sample VCF/BCF files (`output.split_multiallelic`).
- Added information to the output VCF/BCF file headers about the ChoCallate version (`##tool`), multi-allelic variant splitting (`##multialleleSplit`), and monomorphic variant removal (`##noMonomorphic`). The consensus threshold header key was renamed to `##consensusThreshold`.

### Removed

- The NM tag was temporarily removed due to issues during merging. It will be restored in future updates. The consensus threshold remains in the VCF/BCF header.

### Fixed

- Fixed a bug that caused GATK to crash when the system's temporary directory became full.

## [2.0.5] - 2026-05-20

### Fixed

- Fixed typos in the MAPPING_BWA and MAPPING_MINIMAP2 processes that caused them to crash when launched.

### Added

- The MERGE_OUTPUTS process now uses a retry strategy with a decreasing number of processors for each attempt.

## [2.0.4] - 2026-05-09

### Changed

- Refractored the pipeline code and file structure towards nf-core standards.

### Fixed

- Fixed hardcoded paths to input files in the config for test run.

## [2.0.3] - 2026-04-27

### Changed

- Removed `maxForks 1` in several steps, including BAM file filtering, BED file generation with coverage depth, duplicate removal, left-alignment of indels, and BCF to VCF conversion. These were performance bottlenecks for running on many samples.
- Cleanup observer for nf-boost has been changed to `'v2'`.

## [2.0.2] - 2026-04-25

### Added

- Implemented launch in a Docker container.

## [2.0.1] - 2026-04-22

### Fixed

- Name collision was fixed in the coverage generation process, arising from using the same `NO_FILE` for both include_bed and exclude_bed.
- During large file processing, the system temporary directory could overflow, which would result in the pipeline crashing. For each process, the system `$TMPDIR` variable is overridden with a temporary directory that uses the process's `$PWD`.

## [2.0.0] - 2026-04-21

### Changed

- Sequential BCF/VCF indexing was replaced with parallel indexing when merging output files for each sample into a single multi-sample file
- The version number was changed to 2.0.0 due to the removal of backward compatibility in version 1.1.0 compared to 1.0.1, in accordance with semantic versioning

## [1.1.0] - 2026-04-19

### Added

- The `include bed` and `exclude bed` arguments have been added as a more flexible replacement for the `custom bed` argument

### Changed

- All test run data has been moved to a separate `test_run` directory to clean up the project root folder

### Removed

- Argument `custom_bed` has been removed

## [1.0.1] - 2026-04-16

### Fixed

- Added the missing validation of the `reference_index_dir` parameter, which was declared as required, but its absence caused a failure during reads mapping instead of an early termination of the pipeline

## [1.0.0] - 2026-04-14

Initial release
