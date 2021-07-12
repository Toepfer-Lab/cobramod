======
Readme
======

.. image:: img/logo.png
  :width: 600
  :align: center

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

Users will receive proper warnings if any of criteria is violated. In the end,
a log is saved so users can have a overview.

Additionally, CobraMod uses the visualization tool `Escher
<https://escher.readthedocs.io/en/latest/>`_ to illustrate flux distrubutions
for given pathway.

Installation/Development
------------------------

Install CobraMod through pip ::

  pip install cobramod

To download the development branch, clone this repository and install it with
the argument :code:`-e`::

  pip install -e .

Additionally, a conda environment file is supplied in *environment.yml*

Functions
---------

- Retrieve data from a database :func:`cobramod.get_data`
- Transform data into COBRApy objects :func:`cobramod.create_object`
- Add reactions from multiple sources :func:`cobramod.add_reactions.`
- Add metabolites from multiple sources :func:`cobramod.add_metabolites`
- Retrieve complete pathways :func:`cobramod.add_pathway`

Check the documentation of the function for more information.
