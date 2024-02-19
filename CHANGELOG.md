# Changelog

## [Unreleased]

### Planned

- Deprecation of SolCyc

## [1.2.0-dev]

### Added

- Memote example files added to docs
- curation logs are now separated by date

### Changed

- Escher is an optional dependency. Pip installation files reflect this change
- Pre-hooks uses now ruff instead of flake8 and black. This is reflected also
in the dev-environment
- Jupyter files for the documentation reflect current changes
- Dropped support for python 3.8 and 3.9
- Logging format is more human-friendly
- Dropped support of mypy. Black is replaced with ruff
- New pyproject.toml
- Test data updated (19.02.24)

### Removed

- Deprecated config files
- Deprecated guide
- Deprecation of SolCyc. Database is not responding
