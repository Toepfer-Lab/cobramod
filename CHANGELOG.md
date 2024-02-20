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
- Usage of credentials for BioCyc database
- New structure `cobramod.retrieval.Data`. It is in charge of parsing into
  different objects
- Extra functions in utils

### Changed

- Set object bug fixed [https://github.com/Toepfer-Lab/cobramod/issues/50]
- Escher is an optional dependency
- Pre-hooks uses now ruff instead of flake8 and black.
- Behaviour to open files if escher is not installed changed to not open browsers
- Dropped support for python 3.8 and 3.9
- Logging format is more human-friendly
- Dropped support of mypy. Black is replaced with ruff
- New pyproject.toml
- Test data updated (19.02.24)
- Function `available_databases` replaced by `cobramod.retrieval.Databases`
- Structure for retrieval refactored and new parsing sub-packaged created
- DataVersion moved to parsing package
- Updated structure of test module
- Removed multiple "private" methods and functions. 
- Docs files updated to these new changes
- Test updated to uses new changes
- Updating functions to not use deprecated pandas functions

### Removed

- Deprecated config files
- Deprecation warning for SolCyc (Database appears to be shutdown)
