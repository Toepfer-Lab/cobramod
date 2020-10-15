Cobramod
========

Cobramod is a extension package based on `COBRApy <https://github.com/
opencobra/cobrapy>`_. This package facilitates the retrieval of metabolic data
from multiple databases such as `Biocyc <https://biocyc.org/>`_ or `KEGG
<https://www.kegg.jp/>`_.

Cobramod is capable of transforming data from these resources into proper
COBRApy objects, such as Reactions and Metabolites and append them to proper
metabolic models.

..
    comment Insert here a graph

This package takes into consideration the reversibility of the reactions,
coefficient of metabolites, mass balance, among others.

The API of cobramod was created in order to facilitate users extend and test
multiple metabolic pathways, while also maintaining the idea of manual
curation. All changes all recorded in a log file so users can monitor the
changes manually.
