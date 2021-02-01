======
Guides
======

Retrieving data
===============

One basic example is to retrieve data from a list of Biocyc identifiers.
CobraMod utilizes `pathlib's <https://docs.python.org/3/library/
pathlib.html>`_ API to create system paths that can be used in multiple opering
systems. A single loop can manage to retrieve the data.

.. code-block:: python

    from cobramod import get_data
    from pathlib import Path


    dir_data = Path.cwd().joinpath("data")
    identifiers = [
        "CPD-14074",
        "CPD-14075",
        "CPD-14076",
        "CPD-14553",
        "CPD-15317",
        "CPD-15322",
        "CPD-15323",
        "CPD-15326"]

    for single in identifiers:
        get_data(directory=dir_data, identifier=single, database="META")

The first argument of :func`cobramod.get_data` represent the system path where
the data will be stored. In this scenario, it symbolises
``$(CurrentWorkingDirectory)/data``.

The next argument indicates the original identifier to retrieve. Lastly, it is
important to mention the database for these identifiers. A corresponding
directory will be created according the database to store locally the new
files. Following these steps, a new directory will be located in the
following tree::

    data
    `-- META
        |-- CPD-14074.xml
        |-- CPD-14075.xml
        |-- CPD-14076.xml
        |-- CPD-14553.xml
        |-- CPD-15317.xml
        |-- CPD-15322.xml
        |-- CPD-15323.xml
        `-- CPD-15326.xml

As expected, a new directory was created to store locally the data. The name
*META* is derived from `Metacyc <https://metacyc.org/>`_, the accumulation of
Biocyc' sub--databases. To see the working databases, load
:data:`cobramod.available_databases`, e.g:

.. code-block:: python

    from cobramod import available_databases

    >>> print(available_databases)
    ['META', 'ARA', 'KEGG', 'BIGG']

Converting data to objects
==========================

One of the goals of CobraMod is to simplify the creation of COBRApy objects,
such as *Metabolites* and *Reactions*. Single objects can be transform
immediately using the function :func:`cobramod.create_object`:

.. code-block:: python

    from cobramod import create_object
    from pathlib import Path

    dir_data = Path.cwd().joinpath("data")

    new_object = create_object(
        identifier="C00026", directory=dir_data,
        database="KEGG", compartment="c")


    >>> print(type(new_object))
    <class cobra.core.metabolite.Metabolite>

In this example, the KEGG's metabolite `C00026
<https://www.genome.jp/dbget-bin/www_bget?C00026>`_ (2-Oxoglutarate), is
identified as a metabolite and is automatically built as a COBRApy object. It
is not necessary to use *get_data* before, as this function obtains the
data on its own.

Adding Objects
==============

CobraMod uses :mod:`cobra`'s  native objects :class:`cobra.Reaction` and
:class:`cobra.Metabolites`. Furthermore, CobraMod includes an extra class
:class:`cobramod.Pathway`, which inherits and expands the attributes and
methods from his class parent :class:`cobra.core.group.Group`.

Metabolites
"""""""""""

To add metabolites to a Model, you can simply use the function
:func:`cobramod.add_metabolites`.

.. code-block:: python

    from cobramod import add_metabolites
    from cobramod.test import textbook_biocyc
    from pathlib import Path

    dir_data = Path.cwd().joinpath("data")

    test_model = textbook_biocyc.copy()
    add_metabolites(
        model=test_model, obj="MET, c", directory=dir_data, database="META"
    )

    >>> print(type(test_model.metabolites.get_by_id("MET_c")))
    <class cobra.core.metabolite.Metabolite>

In this example, a copy for a test model is made. *add_metabolites* needs as
the first argument the model to extend. The *obj* represent, in this case,
a string with the identifier and the corresponding compartment.

The syntax for a custom metabolite:

:code:`formatted_identifier, name, compartment, chemical_formula,
molecular_charge`

Otherwise, retrieve from database with:

:code:`metabolite_identifier, compartment`

Moreover, instead of a single string you can also add a list with strings. For
the next example, a list is goes to the argument *obj* and includes a new
identifier:

.. code-block:: python

    add_metabolites(
        model=test_model,
        obj=["MET, c", "SUCROSE, c"],
        directory=dir_data,
        database="META",
    )
    >>> print(type(test_model.metabolites.get_by_id("MET_c")))
    <class cobra.core.metabolite.Metabolite>

    >>> print(type(test_model.metabolites.get_by_id("SUCROSE_c")))
    <class cobra.core.metabolite.Metabolite>

There is also the option to instead give the path for a file with the strings.
For instance, given the file *metabolites.txt* in the current working directory
with the content:

.. code-block::

    SUCROSE, c
    MET, c
    MALTOSE_c, MALTOSE[c], c, C12H22O11, 1

You can include the file in the *obj* argument.

.. code-block:: python

    dir_data = Path.cwd().joinpath("data")
    file = Path.cwd().joinpath("metabolites.txt")
    test_model = textbook_biocyc.copy()
    >>> print(len(test_model.metabolites))
    72

    add_metabolites(
        model=test_model, obj=file, directory=dir_data, database="META",
    )
    >>> print(len(test_model.metabolites))
    75

Additionally, regular COBRApy Metabolite can be added (lists are also
supported):

.. code-block:: python

    from cobramod import add_metabolites
    from cobramod.test import textbook, textbook_biocyc

    metabolite = textbook.metabolites.get_by_id("xu5p__D_c")
    test_model = textbook_biocyc.copy()
    add_metabolites(model=test_model, obj=metabolite)

    >>> print(type(test_model.metabolites.get_by_id("xu5p__D_c")))
    <class cobra.core.metabolite.Metabolite>

Reactions
"""""""""

Very much as adding metabolites, a model can be also extended with reactions.
CobraMod includes the function :func:`cobramod.add_reactions`.

.. code-block:: python

    from cobramod.test import textbook_kegg
    from cobramod import add_reactions
    from pathlib import Path

    dir_data = Path.cwd().joinpath("data")
    test_model = textbook_kegg.copy()

    add_reactions(
        model=test_model, obj="R04382, c", directory=dir_data, database="KEGG"
    )

    >>> print(type(test_model.reactions.get_by_id("R04382_c")))
    <class cobra.core.reaction.Reaction>

The first argument represents the model to add the reactions. The next argument
*obj* can pass a string with the identifier and the compartment for the
location to take place. Then, the data directory and the database must be
passed.

Moreover, a list with multiple string can be passed as well.

.. code-block:: python

    add_reactions(
        model=test_model,
        obj=["R04382, c", "R02736 ,c"],
        directory=dir_data,
        database="KEGG",
    )

    >>> print(type(test_model.reactions.get_by_id("R04382_c")))
    <class cobra.core.reaction.Reaction>

Additionally, instead of using a list, the path of a file can be used as well.
Given the file *reactions.txt* in the current working directory with:

.. code-block::

    R04382, c
    R02736, c
    C06118_ce, digalacturonate transport | C06118_c: -1, C06118_e:1


You can pass a :class:`pathlib.Path` in obj for the represantion of the file.

.. code-block:: python

    from cobramod.test import textbook_kegg
    from cobramod import add_reactions
    from pathlib import Path

    dir_data = Path.cwd().joinpath("data")
    test_model = textbook_kegg.copy()
    file = Path.cwd().joinpath("summary.txt")

    >>> print(len(test_model.reactions))
    95

    add_reactions(
        model=test_model, obj=file, directory=dir_data, database="KEGG",
    )
    >>> print(len(test_model.reactions))
    98

The syntax to retrieve reactions:

:code:`original_identifer, compartment`

in case of custom reactions:

:code:`reaction_identifier, reaction_name | metabolite_identifier1:
coefficient,`
:code:`metabolite_identifier2:coefficient, ..., metabolite_identifierX:
coefficient`

Identifiers of metabolites have to end with an underscore and a
compartment:

:code:`OXYGEN-MOLECULE_c: -1`

Finally, regular COBRApy reactions can be added.

.. code-block:: python

    from cobramod.test import textbook_kegg, textbook
    from cobramod import add_reactions
    from pathlib import Path

    dir_data = Path.cwd().joinpath("data")
    test_model = textbook_kegg.copy()
    reaction = textbook.reactions.get_by_id("ACALDt")

    add_reactions(model=test_model, obj=reaction, directory=dir_data)
    >>> print(type(test_model.reactions.get_by_id("ACALDt")))
    <class cobra.core.reaction.Reaction>

.. note:: If the reaction identifies that either the reaction or its
   metabolites are already in the model, under another name, them these
   already-in-model metabolites will be used instead. This is a security
   behaviour to prevent double metabolites.

Pathways
""""""""

CobraMod can add complete pathways into the metabolic models. Using the
function :func:`cobramod.add_pathway`, either a sequence of reaction
identifiers or a pathway identifier can be used as arguments
the model.

.. code-block:: python

    from pathlib import Path
    from cobramod.test import textbook

    dir_data = Path.cwd().joinpath("data")

    >>> print(textbook.optimize().objective_value)
    0.8739215069684307

The original metabolic model `e_coli_core` from COBRApy is loaded under the
name :obj:`cobramod.test.textbook`. It shows an optimation  value of 0.874. For
this example, the identifier `ACETOACETATE-DEG-PWY
<https://biocyc.org/ECOLI/new-image?object=ACETOACETATE-DEG-PWY>`_   will be
used for the test model. This specific pathway has two reactions, in which six
metabolites participates:

.. code-block:: python

    test_model = textbook.copy()
    add_pathway(
        model=test_model,
        pathway="ACETOACETATE-DEG-PWY",
        directory=dir_data,
        database="META",
        compartment="c",
    )
    >>> print(type(test_model.groups.get_by_id("ACETOACETATE-DEG-PWY")))
    <class cobramod.core.pathway.Pathway>

..
    TODO: finish summary and mention the summary in this part.

The pathways included two new metabolites. Thus, sink reactions are
automatically built for them, if needed. In this case, only one sink reaction
is created since the second metabolite can be created from another reaction.

.. code-block:: python

  >>> print(test_model.sinks)
  [<Reaction SK_3_KETOBUTYRATE_c at 0x7f8b1b7bc910>]
  >>> print(test_model.optimize().objective_value)
  20.349250465464955


All the changes, are written into a a record file.

.. code-block::

   ... INFO Data for "ACETOACETATE-DEG-PWY" retrieved.
   ... WARNING Metabolite 'ACETYL-COA' found in given model under 'accoa_c'
   ... WARNING Metabolite 'ACET' found in given model under 'ac_c'
   ... INFO Object 'ACETOACETYL-COA-TRANSFER-RXN' identified as a reaction
   ... WARNING Metabolite 'ACETYL-COA' found in given model under 'accoa_c'
   ... WARNING Metabolite 'CO-A' found in given model under 'coa_c'
   ... INFO Object 'ACETYL-COA-ACETYLTRANSFER-RXN' identified as a reaction
   ... INFO Reaction "ACETOACETYL_COA_TRANSFER_RXN_c" added to model
   ... INFO Testing reaction "ACETOACETYL_COA_TRANSFER_RXN_c"
   ... WARNING Sink reaction created for "3_KETOBUTYRATE_c"
   ... WARNING Sink reaction created for "ACETOACETYL_COA_c"
   ... WARNING Demand reaction for "ACETOACETYL_COA_c" removed
   ... WARNING Demand reaction for "ACETOACETYL_COA_c" removed
   ... INFO Reaction "ACETYL_COA_ACETYLTRANSFER_RXN_c" added to model
   ... INFO Testing reaction "ACETYL_COA_ACETYLTRANSFER_RXN_c"
   ... WARNING Demand reaction for "ACETOACETYL_COA_c" removed
   ... WARNING Sink reaction for "ACETOACETYL_COA_c" removed

It can we be seen that in lines 2,3,5 and 6, CobraMod recognizes metabolites
under another identifiers already in the model, i.e. they will used them
instead.

In this scenario, the objective value changed drastically due to insertion of
a sink reaction. It can be seen that both reaction are being activated if
their fluxes are checked.

.. code-block:: python

  >>> print(f.fluxes["ACETOACETYL_COA_TRANSFER_RXN_c"])
  838.8592591516366
  >>> print(f.fluxes["ACETYL_COA_ACETYLTRANSFER_RXN_c"])
  -838.8592591516366

If the sink reaction *SK_3_KETOBUTYRATE_c* gets removed, the fluxes for this
new pathways are deactivated since there is no reaction to synthetize the
start-metabolite.

.. code-block:: python

  test_model.remove_reactions(["SK_3_KETOBUTYRATE_c"])
  f = test_model.optimize()
  >>> print(f.fluxes["ACETOACETYL_COA_TRANSFER_RXN_c"])
  0.0
  >>> print(f.fluxes["ACETYL_COA_ACETYLTRANSFER_RXN_c"])
  0.0
  >>> print(f.objective_value)
  0.873921506968428

Similar results can be achieved using a sequence. For this example, three
reactions from the `mixed acid fermentation
<https://biocyc.org/META/NEW-IMAGE?type=PATHWAY&object=FERMENTATION-PWY>`_
pathway from Metacyc will be added to the metabolic model:

.. code-block:: python

    from pathlib import Path
    from cobramod import add_pathway
    from cobramod.test import textbook_biocyc

    dir_data = Path.cwd().joinpath("data")
    test_model = textbook_biocyc.copy()
    sequence = ["PEPDEPHOS-RXN", "PYRUVFORMLY-RXN", "FHLMULTI-RXN"]

    >>> print(len(test_model.reactions))
    95

    add_pathway(
        model=test_model,
        pathway=sequence,
        directory=dir_data,
        database="META",
        compartment="c",
        group="test_group"
    )
    >>> print(len(test_model.reactions))
    99

We defined the argument *group* as `test_group`. If we search for that
identifier in the model we will find the new group with these reactions as
members:


.. code-block:: python

    >>> print(type(test_model.groups.get_by_id("test_group")))
    <class cobramod.core.pathway.Pathway>

    >>> print(test_model.groups.get_by_id("test_group").members)
    [<Reaction PEPDEPHOS_RXN_c at 0x7faed876af50>,
    <Reaction PYRUVFORMLY_RXN_c at 0x7faedb8b7c10>,
    <Reaction FHLMULTI_RXN_c at 0x7faed8703a90>]

As expected as the prior example, a extra sink reaction was created since there
is no hydrogen metabolite in the model:

.. code-block:: python

  >>> print(test_model.sinks)
  [<Reaction SK_HYDROGEN_MOLECULE_c at 0x7fb1ff2897d0>]

  >>> print(test_model.optimize().objective_value)
  0.8739215069684305
