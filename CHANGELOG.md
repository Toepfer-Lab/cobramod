# Changelog

## [Unreleased]

### Planned

- Parallelization of downloads and cross-references
- Visualization with other tools than Escher
- Deprecation of SolCyc
- Deprecation of PairDictionary
- Deprecation of Escher's official binding 

## [1.3.0]

### Added

- A new Escher binding independent of the original one [Escher-custom]
- A new Visualization option that displays cobra Objects as a Force Directed Graph
- Introduced [Galata](https://github.com/jupyterlab/jupyterlab/tree/main/galata) as a new test suite
  - Galata automates testing of visualisation as it can compare reference images with outpus created in Jupyter Notebooks

### Changed

- Type checks
  - Reintroduced MyPy as Ruff does not check types according to [ruff's FAQ](https://docs.astral.sh/ruff/faq/#how-does-ruff-compare-to-mypy-or-pyright-or-pyre)
- GitHub Actions
  - GitHub Actions now run tests on Windows and macOS once again 
  - New workflow for the creation of Galata reference images
  - GitHub Action has been divided into smaller sections to provide a better overview of failed steps
- Documentation
  - Now uses [Furo](https://github.com/pradyunsg/furo?tab=readme-ov-file) as a theme
  - Provides additional live examples of the new visualization options
- Tests
  - Tests no longer check for specific database versions, this is only done in a controlled manner for the db_version module
  - Tests for visualization using Escher now use Escher-Custom
  - Several visualization tests were replaced by an equivalent test using Galata
- environment.yml
  - Some dependencies now use their conda-forge package to reduce build time.
  - 'sphinx-autoapi' was set to a dev version to be able to create the documentation using python 3.12
- Maintenance
  - Several adjustments due to announced deprecations

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
