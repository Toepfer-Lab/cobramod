======
Guides
======

Retrieving data
===============

One basic example is to retrieve data for a list of Biocyc identifiers.
Cobramod utilizes `pathlib's <https://docs.python.org/3/library/
pathlib.html>`_ API to create system paths that can be used in multiple opering
systems. A single for loop manage to gather the data from Biocyc.

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

The first argument of *get_data* represent the system path where the data will
be stored. In this scenario, it symbolises
``$(CurrentWorkingDirectory)/data``.

The next argument indicates the original identifier to retrieve. Lastly, it is
important to mention the database for these identifiers. A corresponding
directory will be created according the database to store locally the new
files. Following these steps, a new directory will be located with the
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
Biocyc' sub--databases. Some working examples for the argument `database`:

    * META
    * ARA
    * KEGG

Converting data to objects
==========================

One of the goals of cobramod in to simplify the creation of COBRApy objects.
Single objects can be transform immediately using the method *create_object*:

.. code-block:: python

    from cobramod import create_object
    from pathlib import Path

    dir_data = Path.cwd().joinpath("data")

    new_object = create_object(
        identifier="C00026", directory=dir_data,
        database="KEGG", compartment="c")

.. code-block:: python

    >>> print(type(new_object))
    <class 'cobra.core.metabolite.Metabolite'>

It is not necessary to use *get_data* before. This method is able to obtain the
data and convert it to a proper object individually.

Adding objects from single filename
===================================

Another possibily to add objects into a metabolic model is to use the proper
methods *add_meta_from_file* and *add_reactions_from_file*.

Metabolites
"""""""""""

Thee are two options to create metabolite objects. either from a database or
built customized. For instance, given the following file *new_metabolites.txt*:

.. code-block::
    :emphasize-lines: 3

    SUCROSE, c
    MET, c
    MALTOSE_c, MALTOSE[c], c, C12H22O11, 1

The first two lines represent identifiers from Metacyc, while the highlighted
line is a custom metabolite. For more information about the format, check the
documentation of :func:`cobramod.meta_string_to_model` ::

    from cobra.test import create_test_model
    from cobramod import add_meta_from_file
    from pathlib import Path


    dir_data = Path.cwd().joinpath("data")
    file_path = Path.cwd().joinpath("new_metabolites.txt")

    test_model = create_test_model(model_name="mini")

    >>> print(len(test_model.metabolites))
    23

    add_meta_from_file(
        model=test_model, filename=file_path, directory=dir_data,
        database="META"
    )
    >>> print(len(test_model.metabolites))
    26

The first argument includes the model to be modified. The second and third
argument represent the path of the filename and the directory to store the
data, respectively. Lastly, add the *database* argument to store locally the
new data.

Reactions
"""""""""
As with the previous method, *add_reactions_from_file* is able to add custom
reactions or build reactions from a database. Given the following file
*new_reactions.txt*:

.. code-block::
    :emphasize-lines: 3

    R04382, c
    R02736, c
    C06118_ce, digalacturonate transport | C06118_c: -1, C06118_e:1

The first lines represent reactions identifiers for the databse KEGG. It is
possible to reaction custom reactions included in the file. For more
information, check the documentation of
:func:`cobramod.add_reactions_from_file`::

    from cobra.test import create_test_model
    from cobramod import add_reactions_from_file
    from pathlib import Path


    dir_data = Path.cwd().joinpath("data")
    file_path = Path.cwd().joinpath("new_reactions.txt")

    test_model = create_test_model(model_name="mini")

    >>> print(len(test_model.reactions))
    18

    add_reactions_from_file(
        model=test_model, filename=file_path, directory=dir_data,
        database="KEGG"
    )
    >>> print(len(test_model.reactions))
    21

Similar to *add_meta_from_file*, the only arguments needed are the metabolic
model to modify, the file path of the reactions, the directory to store the
data and the name of the database.

Adding Pathways
"""""""""""""""

.. note::
    Currently, only the pathway syntax of Metacyc is working.

Cobramod can add complete pathways into the metabolic models. Using the method
:func:`cobramod.add_graph_to_model`, either a sequence of reaction identifiers
or the original pathway identifier for a database can be used to be added into
the model.

.. code::

  from pathlib import Path
  from cobramod import add_graph_to_model
  from cobramod.test import mini_model

  dir_data = Path.cwd().joinpath("data")

  >>> print(mini_model.optimize().objective_value)
  0.8739215069684307

The original metabolic model `e_coli_core` from COBRApy shows an optimation
value of 0.874. For this example, the identifier `ACETOACETATE-DEG-PWY
<https://biocyc.org/ECOLI/new-image?object=ACETOACETATE-DEG-PWY>`_   will be
used for the test model. This specific pathway has two reactions, in which six
metabolites participates.::

  test_model = mini_model.copy()
  >>> add_graph_to_model(
         model=test_model,
         graph="ACETOACETATE-DEG-PWY",
         directory=dir_data,
         database="META",
         compartment="c",
      )
  --------------------
  Model: e_coli_core
  Original attributes:
  Reactions: 95
  Metabolites: 72
  Boundary reactions 20
  --------------------
  New attributes:
  Reactions: 98
  Metabolites: 74
  Boundary reactions: 21
  --------------------

The output of the method is a short summary about the change of attributes for
the model. The pathways included two metabolites, which were not in the model
and thus, sink reactions are automatically built for them. However, only one
sink reaction is created since the second metabolite can be created from
another reaction. As expected. a total of three new reactions are added, from
which one is a sink::

  >>> print(test_model.sinks)
  [<Reaction SK_3_KETOBUTYRATE_c at 0x7f8b1b7bc910>]
  >>> print(test_model.optimize().objective_value)
  20.349250465464955


All the changes, are written into a a record file.

.. code-block:: text


  2020-10-28 15:18:30,164 INFO Data for "ACETOACETATE-DEG-PWY" retrieved
  2020-10-28 15:18:30,168 INFO Data for "ACETOACETYL-COA-TRANSFER-RXN"\
  retrieved.
  2020-10-28 15:18:30,181 INFO Data for "ACETYL-COA-ACETYLTRANSFER-RXN"
  retrieved.
  2020-10-28 15:18:30,208 INFO Reaction "ACETOACETYL_COA_TRANSFER_RXN_c" added
  to model
  2020-10-28 15:18:30,208 INFO Testing reaction
  "ACETOACETYL_COA_TRANSFER_RXN_c"
  2020-10-28 15:18:30,214 WARNING Sink reaction created for "3_KETOBUTYRATE_c"
  2020-10-28 15:18:30,216 WARNING Sink reaction created for "ACETOACETYL_COA_c"
  2020-10-28 15:18:30,217 WARNING Demand reaction for "ACETOACETYL_COA_c"
  removed
  2020-10-28 15:18:30,219 WARNING Demand reaction for "ACETOACETYL_COA_c"
  removed
  2020-10-28 15:18:30,222 INFO Reaction "ACETYL_COA_ACETYLTRANSFER_RXN_c"
  added to model
  2020-10-28 15:18:30,222 INFO Testing reaction
  "ACETYL_COA_ACETYLTRANSFER_RXN_c"
  2020-10-28 15:18:30,225 WARNING Demand reaction for "ACETOACETYL_COA_c"
  removed
  2020-10-28 15:18:30,226 WARNING Sink reaction for "ACETOACETYL_COA_c" removed

In this scenario, the objective value changed drastically due to insertion of
the sink reaction. It can be seen that both reaction are being activated if
their fluxes are checked::

  >>> print(f.fluxes["ACETOACETYL_COA_TRANSFER_RXN_c"])
  838.8592591516366
  >>> print(f.fluxes["ACETYL_COA_ACETYLTRANSFER_RXN_c"])
  -838.8592591516366

If the sink reaction, in this case *SK_3_KETOBUTYRATE_c* gets removed, the
fluxes for this new pathways are deactivated since there is no reaction to
synthetize the start-metabolite::

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
pathway from Metacyc will be added to the metabolic model::


  from pathlib import Path
  from cobramod import add_graph_to_model
  from cobramod.test import mini_model

  dir_data = Path.cwd().joinpath("data")
  test_model = mini_model.copy()

  sequence = ["PEPDEPHOS-RXN", "PYRUVFORMLY-RXN", "FHLMULTI-RXN"]
  >>> add_graph_to_model(
          model=test_model,
          graph=sequence,
          directory=dir_data,
          database="META",
          compartment="c",
      )
  Model: e_coli_core
  Original attributes:
  Reactions: 95
  Metabolites: 72
  Boundary reactions 20
  --------------------
  New attributes:
  Reactions: 99
  Metabolites: 74
  Boundary reactions: 21
  --------------------

As expected as the prior example, a extra sink reaction was created since there
is no hydrogen metabolite in the model::

  >>> print(test_model.sinks)
  [<Reaction SK_HYDROGEN_MOLECULE_c at 0x7fb1ff2897d0>]

  >>> print(test_model.optimize().objective_value)
  0.8739215069684305
