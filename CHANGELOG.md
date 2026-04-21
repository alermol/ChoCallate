# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
