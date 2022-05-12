![Generic badge](https://img.shields.io/badge/python-3.7%20%7C%203.8%20%7C%203.9%20%7C%203.10-blue)
![GitHub](https://img.shields.io/github/license/Toepfer-Lab/cobramod)
![Read the Docs (version)](https://img.shields.io/readthedocs/cobramod/latest)
![Tests](https://img.shields.io/github/workflow/status/Toepfer-Lab/cobramod/Test%20build%20and%20publish%20Cobramod%20to%20PyPI?label=tests)
![Version](https://img.shields.io/pypi/v/cobramod?label=version)
[![Downloads](https://pepy.tech/badge/cobramod)](https://pepy.tech/project/cobramod)


CobraMod: A pathway-centric curation tool for constraint-based metabolic models
===============================================================================

![image](https://raw.githubusercontent.com/Toepfer-Lab/cobramod/master/docs/source/img/logo.png)

CobraMod is a Python 3 open-source package for pathway-centric curation
of genome-scale metabolic models (GEMs). It builds upon the [COBRApy
toolbox](https://opencobra.github.io/cobrapy/) and offers a
comprehensible set of functions for semi-automated network extension,
curation and visualization. CobraMod supports all databases from the
[BioCyc collection](https://biocyc.org/), the [KEGG
database](https://www.genome.jp/kegg/), and the [BiGG Models
repository](http://bigg.ucsd.edu/) and can directly interact with Escher
for pathway and flux visualization.

This package converts pathway information into native COBRApy objects
and quality-checks them while adding them to the model. This includes
testing for:

-   duplicate elements
-   correct chemical formula
-   mass balance of reactions
-   reaction reversibility
-   capability to carry non-zero fluxes
-   adding available gene information
-   MEMOTE compliance
-   available cross-references

CobraMod offers user-friendly tracking of the curation process with
summary output and log files and customized pathway and flux
visualization with [Escher](https://escher.github.io/).

Installation
------------

CobraMod can easily be installed using pip :

    pip install cobramod

Functions
---------

This package offers multiple functions for modifying and extending GEMs:

-   Retrieve metabolic pathway information from a database
    `cobramod.get_data`
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
    author = {Camborda, Stefano and Weder, Jan-Niklas and Töpfer, Nadine},
    title = "{CobraMod: a pathway-centric curation tool for constraint-based metabolic models}",
    journal = {Bioinformatics},
    volume = {38},
    number = {9},
    pages = {2654-2656},
    year = {2022},
    month = {02},
    abstract = "{COnstraint-Based Reconstruction and Analysis of genome-scale metabolic models has become a widely used tool to understand metabolic network behavior at a large scale. However, existing reconstruction tools lack functionalities to address modellers' common objective to study metabolic networks on the pathway level. Thus, we developed CobraMod—a Python package for pathway-centric modification and extension of genome-scale metabolic networks. CobraMod can integrate data from various metabolic pathway databases as well as user-curated information. Our tool tests newly added metabolites, reactions and pathways against multiple curation criteria, suggests manual curation steps and provides the user with records of changes to ensure high quality metabolic reconstructions. CobraMod uses the visualization tool Escher for pathway representation and offers simple customization options for comparison of pathways and flux distributions. Our package enables coherent and reproducible workflows as it can be seamlessly integrated with COBRApy and Escher.The source code can be found at https://github.com/Toepfer-Lab/cobramod/ and can be installed with pip. The documentation including tutorials is available at https://cobramod.readthedocs.io/.}",
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

You can contribute to CobraMod by cloning the repository and installing
it in developer mode using pip :

    pip install -e .

A conda environment file is supplied (*environment.yml*). This file
contains all the dependencies to ensure reproducibility of the package.
We encourage pull requests. CobraMod uses unit testing and new tests are
welcome.

To report bugs and suggestions, please create an issue using the
corresponding tags at <https://github.com/Toepfer-Lab/cobramod/issues>.
