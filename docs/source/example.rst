======
Guides
======

Retrieving data
===============

One basic example is to retrieve data for a list of Biocyc identifiers.
Cobramiod utilizes `pathlib's <https://docs.python.org/3/library/
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

It is possible to build metabolites either from a database or customized.
For instance, given the following file *new_metabolites.txt*:

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
