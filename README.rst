Cobram d
========
.. image:: docs/source/img/logo.png
  :width: 400
  :align: center

Cobramod is a extension package based on `COBRApy <https://github.com/
opencobra/cobrapy>`_. This package facilitates the retrieval of metabolic data
from multiple databases such as `Biocyc <https://biocyc.org/>`_ or `KEGG
<https://www.kegg.jp/>`_.

actuate
Cobramod is capable of transforming data from these resources into proper
COBRApy objects, such as reactions and metabolites while appending them to
proper metabolic
models.

This package takes into consideration the reversibility of the reactions,
coefficient of metabolites, mass balance, among others attributes.

The API of cobramod was created to help users extend and test multiple
metabolic pathways, while maintaining the idea of manual curation. All changes
all recorded in a log file that users can monitor the changes manually.
