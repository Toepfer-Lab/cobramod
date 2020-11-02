Cobramod
========
.. image:: docs/source/img/logo.png
  :width: 400
  :align: center

Cobramod is a python extension package based on `COBRApy <https://github.com/
opencobra/cobrapy>`_. This package facilitates the retrieval of metabolic data
from multiple databases such as `Biocyc <https://biocyc.org/>`_ or `KEGG
<https://www.kegg.jp/>`_.

Cobramod is capable of transforming data from these resources into proper
COBRApy objects, such as reactions and metabolites while appending them to
proper metabolic models, such that users are able to track changes in the
metabolic model.

This package takes into consideration the reversibility of the reactions,
coefficient of metabolites, mass balance, among others attributes.

The API of cobramod was build to help users extend and test multiple
metabolic pathways, while maintaining the idea of manual curation. All changes
all recorded in a log file, which users can monitor changes manually.
