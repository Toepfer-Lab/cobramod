dbwalker
*************************

.. figure:: /img/dbwalker.svg
   :scale: 50 %
   :alt: The workflow for acquiring additional identifiers



.. panels::

    Starting with CobraMod 2.0, the addition of extra identifiers from databases that aren't the original source has been significantly improved. Instead of relying on MetaNetX in general, identifiers from other databases are now identified using InChI, InChIKey, and SMILES. To do this, they are first queried in the original database and used to identify equivalent database-specific identifiers in other databases. Since errors sometimes occur in other tools, depending on how accurately they perform the matching, CobraMod includes a verification process. On the one hand, within a database, database-unspecific identifiers must refer to database-specific identifiers. On the other hand, all database-unspecific identifiers must point to the same database-specific identifier.

    ---

    .. csv-table:: Databases directly supported for identifier queuing
      :header: "Database", "Metabolites", "Reactions"
      :width: 65%
      :widths: auto
      :align: center

      "BioCyc", "x", ""
      "PubChem", "x", ""
      "KEGG", "x", ""
      "ChEBI", "x", ""

      "BIGG", "", ""
      "enviPath",  "", ""
      "HMDB",  "", ""
      "LipidMaps",  "", ""
      "Reactome",  "", ""
      "Rhea",  "", ""
      "SABIO-RK",  "", ""
      "SwissLipids", "", ""
      "The SEED",  "", ""

----------
BIGG
----------

Bigg does not provide independent mappings for entries within their database, instead they incorporate information
originally provided by MetaNetX. As such identifiers already in present for example pointing towards BioCyc are already
handled by the BioCyc implementation. Mappings to Bigg an from Bigg to other databases will either way go be constructed
by MetaNetX so we dont provide a specific implementation for Bigg.
