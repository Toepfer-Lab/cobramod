==============================================================================
CobraMod: A pathway-centric curation tool for contraint-based metabolic models
==============================================================================

CobraMod is an open-source Python 3 package intended to be an extension of
`COBRApy <https://github.com/opencobra/cobrapy>`_. This package facilitates
the retrieval of biochemical
data from multiple databases such as `Biocyc <https://biocyc.org/>`_, `KEGG
<https://www.kegg.jp/>`_ and `BiGG <http://bigg.ucsd.edu/>`_; and enables
extending metabolic models.

It is capable of transforming data from these resources into
proper COBRApy objects, such as :class:`cobra.Reaction` and
:class:`cobra.Metabolite` while appending
them singlely or as sets, denominated `Pathways`. CobraMod focus on the
curation of these sets and takes into consideration following multiple
criteria:

- chemical formula
- duplicate elements
- reversibility and mass balance of reactions
- capability to carry non‑zero fluxes
- changes in the model.

Users will receive proper warnings if any of criteria are violated. In the end,
a log is saved so users can have a overview.

Additionally, CobraMod uses the visualization tool `Escher
<https://escher.readthedocs.io/en/latest/>`_ to illustrate flux distrubutions
for given pathway.

.. toctree::
   :maxdepth: 2
   :numbered:
   :caption: Contents

   README.rst
   functional.rst
   installation.rst
   how_to.ipynb
   module/index.html

.. toctree::
   :maxdepth: 0
   :caption: Test case
   :numbered:

   Glucosinolate scenario <GLS.ipynb>

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
