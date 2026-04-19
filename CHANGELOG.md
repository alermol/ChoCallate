# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
