# Changelog

## [Unreleased]

### Planned

- Parallelization of downloads and cross-references
- Visualization with other tools than Escher
- Deprecation of SolCyc
- Deprecation of PairDictionary

## [1.2.0]

### Added

- Memote example files added to docs
- Curation and visualization logs are now separated by date in a 'logs' dir
- Usage of credentials for Biocyc database

### Changed

- Escher is an optional dependency. Pip installation files reflect this change
- Pre-hooks uses now ruff instead of flake8 and black. This is reflected also
in the dev-environment
- Jupyter files for the documentation reflect current changes
- Behaviour to open files if escher is not installed changed to not open browsers
- Dropped support for python 3.8 and 3.9
- Logging format is more human-friendly
- Dropped support of mypy. Black is replaced with ruff
- New pyproject.toml
- Test data updated (19.02.24)

### Removed

- Deprecated config files
- Deprecated guide
- Deprecation of SolCyc. Database appears to be shutdown
