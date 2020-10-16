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

One of the goals of cobramod in to simplify the creation ob COBRApy objects.
Single objects can transform immediately using the method *create_object*:

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
