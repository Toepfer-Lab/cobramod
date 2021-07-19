.. image:: https://img.shields.io/github/license/Toepfer-Lab/cobramod
   :alt: GitHub
.. image:: https://img.shields.io/readthedocs/cobramod/latest
   :alt: Read the Docs (version)

===============================================================================
CobraMod: A pathway-centric curation tool for constraint-based metabolic models
===============================================================================

.. image:: docs/source/img/logo.png
  :width: 600
  :align: center

CobraMod is a Python 3 open-source package for pathway-centric curation of
genome-scale metabolic models (GEMs). It builds upon the
`COBRApy toolbox <https://opencobra.github.io/cobrapy/>`_
and offers a comprehensible set of functions for semi-automated network
extension, curation and visualization. CobraMod supports all databases from the
`BioCyc collection <https://biocyc.org/>`_, the
`KEGG database <https://www.genome.jp/kegg/>`_, and the
`BiGG Models repository <http://bigg.ucsd.edu/>`_. and can directly interact
with Escher for pathway and flux
visualization.

This package converts pathway information into native COBRApy objects and
quality-checks them before adding them to the model. CobraMod will use and
parse the exact information from the metabolic pathway information. This
includes  testing for:

- duplicate elements
- correct chemical formula according to the data
- assignment of genes
- mass balance of reactions
- reaction reversibility
- capability to carry non-zero fluxes

CobraMod offers user-friendly tracking of the curation process with summary
output and log files and customized pathway and flux visualization with Escher.
CobraMod uses `Escher <https://escher.github.io/>`_ for visualizing pathways
and flux distributions and offers several customization options.

Installation
---------------

CobraMod can easily be installed using pip. ::

  pip install cobramod


Functions
-------------

This package offers multiple functions for modifying and extending GEMs:

- Retrieve metabolic pathway information from a database
  :func:`cobramod.get_data`
- Transform stored data into COBRApy objects :func:`cobramod.create_object`
- Add metabolites from multiple sources :func:`cobramod.add_metabolites`
- Add reactions from multiple sources :func:`cobramod.add_reactions`
- Test reaction capability to carry a non-zero flux
  :func:`cobramod.test_non_zero_flux`
- Add pathway to a model :func:`cobramod.add_pathway`

Check the `documentation <https://cobramod.readthedocs.io/>`_ for more
information.

License
------------
CobraMod is licensed under the GPL-3 License. Read :doc:`LICENSE` for more
information.


Development
-------------------

You can contribute to CobraMod by cloning the repository and installing it in
developer mode using pip::

  pip install -e .

A conda environment file is supplied (*environment.yml*). This file supplies
the environment that  we use to ensure the reproducibility of the package. To
report bugs and suggestions, please create an issue using the corresponding
tags at https://github.com/Toepfer-Lab/cobramod/issues.

We encourage pull requests. CobraMod uses unit testing and new tests are
welcome.
