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

This package converts pathway information into native COBRApy
objects and quality-checks them before adding them to the model. This includes
testing for:

- duplicate elements
- correct chemical formula
- mass balance of reactions
- reaction reversibility
- capability to carry non-zero fluxes
- adding available gene information

CobraMod offers user-friendly tracking of the curation process with summary
output and log files and customized pathway and flux visualization with Escher.
CobraMod uses `Escher <https://escher.github.io/>`_ for visualizing pathways
and flux distributions and offers several customization options.

.. raw:: html

   <table style="width: 100%; border-collapse: collapse; border-style: none;" border="0">
   <tbody>
   <tr>
   <td style="width: 100%; text-align: center;">
   <h4>Example of the visualization method</h4>
   </td>
   </tr>
   <tr>
   <td style="width: 100%; text-align: justify;"><img src="_static/pathway01.png" alt="" /><span style="color: #999999;">This pathway representation uses Escher. We created three fictive reactions and added them to a toy model. The arrows represent the flux of the reactions. Using CobraMod we visualize the pathway and set colors for the fluxes. Red color represents negative fluxes and green positive fluxes.</span></td>
   </tr>
   </tbody>
   </table>

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

- Using COBRApy objects
- Using a text file with database-specific identifiers for the metabolites,
  reactions, or pathways.
- Using a single string with the database-specific identifier of the object
- Using a text file or string with user-defined metabolite,
  reaction information

CobraMod includes the class :class:`cobramod.Pathway` which inherits all
methods of the COBRApy :class:`cobra.core.group.Group` and additionally
includes the method :func:`cobramod.Pathway.visualize` for creating pathway and
flux maps.


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
   test_cases.ipynb
   API <module/cobramod/index>
