Settings
*************************

Cobramod's settings object is the central point for controlling its behaviour, setting local save paths, and managing login information.

--------------------
BioCyc-Login
--------------------

To access BioCyc, a license is required. CobraMod offers three ways to provide the necessary login information: a file containing the username and password, environment variables, or passing it directly to the settings object.

File
--------------------
When using a BioCyc login file, it must have two lines. The first line is the username, and the second line is the password. No other text is allowed.

.. code-block:: text

   <Username>
   <Password>

This file can then be used via a settings object, as shown in the following:

.. code-block:: python

   from cobramod import Settings

   settings = Settings()

   settings.SetBioCycLogin("./fileWithLoginInfo.txt")


.. note::
Cobramod also automatically checks the local folder for a file called 'credentials.txt' and, if found, loads the login information from it.

Environment Variables
---------------------

Instead of using a file, CobraMod also supports environment variables. For this, the username needs to be assigned to a variable called 'BIOCYC_NAME', and the password to a variable called 'BIOCYC_PASSWORD'. If no login information is provided, CobraMod checks automatically for the presence of those variables and uses them when necessary.

.. code-block:: text
    BIOCYC_NAME
    BIOCYC_PASSWORD

Via a Settings object
-----------------------

Last but not least, the login information can also be directly passed to the settings object. This can be done via assigning the login information to the ‘biocyc_name’ and ‘biocyc_password’ variables.

.. code-block:: python

   from cobramod import Settings

   settings = Settings()

   settings.biocyc_name = "someone"
   settings.biocyc_password = "APassword"

--------------------
SMILES as identifier
--------------------
Since SMILES are not recognized by 'identifiers.org' as allowed cross-references, we do not add them to CobraPy objects (Metabolites) by default. We do allow them to be introduced as identifiers, but this behaviour must be explicitly enabled by setting the 'add_smiles_as_cross_reference' boolean in the settings object to True.

.. code-block:: python

   from cobramod import Settings

   settings = Settings()

   settings.add_smiles_as_cross_reference = True

.. note::
   This setting only affects cross-references that are retrieved from MetaNetX after it is changed. Cross-references that are already cached in the data directory (folder 'XRef') are used as they are, with or without SMILES.

-----------------------------
MetaNetX specific prefixes
-----------------------------
MetaNetX returns most identifiers two times: once with the prefix registered in 'identifiers.org' and once with a MetaNetX specific prefix, e.g. 'kegg.compound:C00031' and 'keggC:C00031'. We do not add cross-references with a MetaNetX specific prefix, such as 'keggC' or 'biggM', to CobraPy objects by default. To add them anyway, set the 'add_metanetx_specific_prefixes' boolean in the settings object to True.

.. code-block:: python

   from cobramod import Settings

   settings = Settings()

   settings.add_metanetx_specific_prefixes = True

Example
-------

For D-glucose with the annotation ``{"kegg.compound": "C00031"}``, MetaNetX returns 36 cross-references. The table shows how CobraMod adds the ones for BiGG, KEGG and SEED to the annotation:

.. list-table::
   :header-rows: 1

   * - Returned by MetaNetX
     - Default (False)
     - True
   * - bigg.metabolite:glc__D
     - added as ``bigg.metabolite``
     - added as ``bigg.metabolite``
   * - biggM:glc__D, biggM:M_glc__D, biggM:glc, biggM:glc-D
     - ignored
     - added as ``biggm``
   * - kegg.compound:C00031
     - added as ``kegg.compound``
     - added as ``kegg.compound``
   * - keggC:C00031, keggC:M_C00031
     - ignored
     - added as ``keggc``
   * - seed.compound:cpd00027
     - added as ``seed.compound``
     - added as ``seed.compound``
   * - seedM:cpd00027, seedM:M_cpd00027
     - ignored
     - added as ``seedm``

The other MetaNetX specific prefixes (e.g. 'keggD', 'metacycM', 'sabiorkM' or 'vmhM') are handled the same way. Their values mostly repeat the 'identifiers.org' identifiers, partly with an additional 'M\_' prefix. A few, such as the old BiGG identifiers 'glc' and 'glc-D' above, are only available with the MetaNetX specific prefixes.
