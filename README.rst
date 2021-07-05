.. image:: docs/source/img/logo.png
  :width: 400
  :align: center

CobraMod is an open-source Python 3 package based on `COBRApy <https://github.com/
opencobra/cobrapy>`_. This package facilitates the retrieval of biochemical
data from multiple databases such as `Biocyc <https://biocyc.org/>`_, `KEGG
<https://www.kegg.jp/>`_ and `BiGG <http://bigg.ucsd.edu/>`_ and enables
extending metabolic models.

It is capable of transforming data from these resources into
proper COBRApy objects, such as `Reactions` and `Metabolites` while appending
them singlely or as sets, denominated `Pathways`. CobraMod focus on the
curation of these sets and takes into consideration following multiple
criteria:

- chemical formula
- duplicate elements
- reversibility and mass balance of reactions
- capability to carry non‑zero fluxes
- changes in the model.

In case that any of the criteria is not met or exceptions are found, the users
will receive proper warnings that are showed in the console and stored in a log
file for a proper overview.

Additionally
The API of cobramod was build to help users extend and test multiple
metabolic pathways, while maintaining the idea of manual curation. All changes
all recorded in a log file, which users can monitor changes manually.
