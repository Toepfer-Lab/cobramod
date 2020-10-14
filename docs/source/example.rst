Guides
======

One basic example is to retrieve data for a list of Biocyc identifiers. Cobramod utilizes
`pathlib's <https://docs.python.org/3/library/pathlib.html>`_ API to create system paths
that can be used in multiple opering systems. A single for loop manage to gather the
data from Biocyc. 


.. code-block:: python

    from cobramod.mod_parser import get_data
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
        get_data(
            directory=dir_data,
            identifier=single,
            database="META")


The first argument of `get_data` represent the system path where the data will be stored.
In this scenario, it symbolises :code:`$(CurrentWorkingDirectory)/data`.
The next argument indicates the original identifier to retrieve. Lastly, it is import to 
mention the database for these identifiers. A corresponding directory will be created
according the database to store locally the new files. Following these steps, a new directory will be located with the following tree:

-- META
    |-- CPD-14074.xml
    |-- CPD-14075.xml
    |-- CPD-14076.xml
    |-- CPD-14553.xml
    |-- CPD-15317.xml
    |-- CPD-15322.xml
    |-- CPD-15323.xml
    |-- CPD-15326.xml
 