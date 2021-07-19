============
Introduction
============

CobraMod is a Python 3 open-source package for pathway-centric curation of
genome-scale metabolic models (GEMs). It builds upon the
`COBRApy toolbox <https://opencobra.github.io/cobrapy/>`_
and offers a comprehensible set of functions for semi-automated network
extension, curation and visualization. CobraMod supports all databases from the
`BioCyc collection <https://biocyc.org/>`_, the
`KEGG database <https://www.genome.jp/kegg/>`_, and the
`BiGG Models repository <http://bigg.ucsd.edu/>`_ and can directly interact
with Escher for pathway and flux
visualization.

CobraMod will use and parse the exact information from the metabolic pathway
information. This package converts pathway information into native COBRApy
objects and quality-checks them before adding them to the model. This includes
testing for:

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


.. image:: img/pathway01.png
  :width: 600
  :align: center


This package offers multiple functions for modifying and extending GEMs:

- Retrieve metabolic pathway information from a database
  :func:`cobramod.get_data`
- Transform stored data into COBRApy objects :func:`cobramod.create_object`
- Add metabolites from multiple sources :func:`cobramod.add_metabolites`
- Add reactions from multiple sources :func:`cobramod.add_reactions`
- Test reaction capability to carry a non-zero flux
  :func:`cobramod.test_non_zero_flux`
- Add pathway to a model :func:`cobramod.add_pathway`

Users can add the biochemical data through different methods:

- Using regular COBRApy objects
- Using a text file with the database-specific identifiers for the reactions,
  metabolic or pathways.
- Using a single string with the database-specific identifier of the object

In addition to using databases, the users can include user-curated reactions
and metabolites. These objects are plain text with simple syntax that can be
include in the text files or directly in the functions.

CobraMod includes a new :class:`cobramod.Pathway`. This class inherits
the methods of the original COBRApy :class:`cobra.core.group.Group`.
This pathway-object includes the method
:func:`cobramod.Pathway.visualize` for creating  pathway maps.


For more information see the docstrings of the respective functions using
`help()` or read the :doc:`documentation <module/cobramod/index>`.

To see all metabolic pathway databases that are currently supported by
CobraMod, load and print :obj:`cobramod.available_databases`.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Table of content

   how_to.ipynb
   installation.rst
   shikimate.ipynb
   GLS.ipynb
   API <module/cobramod/index>
