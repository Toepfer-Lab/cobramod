<div align="center">

![Static Badge](https://img.shields.io/badge/python-3.10%7C3.11%7C3.12-%20blue)
![GitHub](https://img.shields.io/github/license/Toepfer-Lab/cobramod)
[![Downloads](https://img.shields.io/pepy/dt/cobramod
)](https://pepy.tech/project/cobramod)

![Version](https://img.shields.io/pypi/v/cobramod?label=version)
![Read the Docs (version)](https://img.shields.io/readthedocs/cobramod/latest)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/Toepfer-lab/cobramod/test-build-and-publish.yml)
![Coverage Status](./docs/source/img/coverage.svg)

CobraMod: A pathway-centric curation tool for constraint-based metabolic models
===============================================================================


![image](https://raw.githubusercontent.com/Toepfer-Lab/cobramod/master/docs/source/img/logo.png)
</div>

CobraMod is a Python 3 open-source package for pathway-centric curation of
genome-scale metabolic models (GEMs). It builds upon the 
[COBRApy toolbox](https://opencobra.github.io/cobrapy/) and offers a 
comprehensible set of functions for semi-automated network extension, curation and visualization.
CobraMod supports all databases from the [BioCyc collection](https://biocyc.org/), 
the [KEGG database](https://www.genome.jp/kegg/), and the [BiGG
Models repository](http://bigg.ucsd.edu/). It optionally can interact with
Escher for pathway and flux visualization.

This package converts pathway information into native COBRApy objects and
quality-checks them while adding them to the model. This includes testing for:

-   duplicate elements
-   correct chemical formula
-   mass balance of reactions
-   reaction reversibility
-   capability to carry non-zero fluxes
-   adding available gene information
-   MEMOTE compliance
-   available cross-references

CobraMod offers user-friendly tracking of the curation process with summary
output with log files and customized pathway and flux visualization with
[Escher](https://escher.github.io/).

Installation
------------

CobraMod can easily be installed using pip :

    pip install cobramod

Functions
---------

This package offers multiple functions for modifying and extending GEMs:

-   Retrieve metabolic pathway information from a database `cobramod.get_data`
-   Transform stored data into COBRApy objects `cobramod.create_object`
-   Add metabolites from multiple sources `cobramod.add_metabolites`
-   Add reactions from multiple sources `cobramod.add_reactions`
-   Test reaction capability to carry a non-zero flux
    `cobramod.test_non_zero_flux`
-   Add pathway to a model `cobramod.add_pathway`
-   Automatic cross-references `cobramod.add_crossreferences`
-   Testing for MEMOTE compliance

Check the [documentation](https://cobramod.readthedocs.io/) for more
information.

Citing CobraMod
---------------

To cite CobraMod, please use the following paper:

[CobraMod: a pathway-centric curation tool for constraint-based metabolic models](https://doi.org/10.1093/bioinformatics/btac119)

```bibtex
@article{10.1093/bioinformatics/btac119,
    author = {Camborda La Cruz, Stefano and Weder, Jan-Niklas and Töpfer, Nadine},
    title = "{CobraMod: a pathway-centric curation tool for constraint-based metabolic models}",
    journal = {Bioinformatics},
    volume = {38},
    number = {9},
    pages = {2654-2656},
    year = {2022},
    month = {02},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btac119},
    url = {https://doi.org/10.1093/bioinformatics/btac119},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/38/9/2654/43481008/btac119.pdf},
}
```
License
-------

CobraMod is licensed under the GPL-3 License. Read the
[LICENSE](https://github.com/Toepfer-Lab/cobramod/blob/master/LICENSE)
for more information.


Development
-----------
CobraMod consists of a Python and JavaScript/TypeScript part. 
The following briefly describes all the commands to build a local 
version from the source files.

### JavaScript

Since we need Node modules to build the Jupyter integrations con CobraMod, 
a package.json is included in the project. All dependencies contained 
in it can be made available locally using [yarn](https://yarnpkg.com/). The foiling command can be used for this:

    yarn install

We use [Vite](https://vitejs.dev/) to build the Javascript part of the Jupyter integrations. So, we can use the following command to build a bundled version of the integrations.

    yarn run vite build

This creates bundled JavaScript files under '/src/cobramod/static'.

### Python

You can contribute to CobraMod by cloning the repository and installing it in
developer mode and using the `dev` dependency group via pip:

    pip install -e ".[dev]"

Optionally, a conda environment file is supplied (*environment.yml*). This
file contains all the dependencies to ensure reproducibility of the package. We
encourage pull requests!

To report bugs and suggestions, please create an issue using the
corresponding tags at <https://github.com/Toepfer-Lab/cobramod/issues>.
