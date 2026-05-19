# Changelog

All notable changes to RDPOWER will be recorded here. This changelog starts
with the May 2026 modernization baseline.

## [R 3.0 / Python 2.0.0] - 2026-05-15

Modernization release prepared across May 15-16, 2026.

### Added

- Added fixed-seed numerical contract coverage and baseline comparison tooling
  for the R, Python, and Stata implementations.
- Added repository release hardening, including GitHub Actions CI, issue and
  pull request templates, Dependabot and security configuration, Git attributes
  and ignore rules, and PyPI trusted publishing support.
- Added modern Python package metadata under `Python/rdpower`, including
  `pyproject.toml`, package-level licensing, README content, and shared utility
  helpers.
- Added scripts for repository validation and Stata help PDF generation.
- Added R summary methods for `rdpower`, `rdmde`, and `rdsampsi` output objects.

### Changed

- Updated the R package to version 3.0, the Python package to version 2.0.0, and
  the Stata distribution date to 2026-05-15.
- Modernized package metadata with current author emails, GPL-3.0 licensing
  metadata, and repository issue links.
- Moved the Python package from `python/rdpower` to `Python/rdpower` and removed
  stale generated Python build artifacts from version control.
- Refreshed README content, R documentation, Python documentation, Stata help
  files and PDFs, and the R, Python, and Stata illustration scripts.
- Aligned R output handling with the new summary methods while preserving the
  package's returned result objects.
- Removed temporary audit baseline files from tracked content after release
  preparation.

### Validated

- Added package-level checks intended to cover R, Python, and Stata numerical
  baselines, repository metadata, Python package tests, and Stata help artifact
  consistency.
