Requirements
============

To run CobraMod, *Python 3.7.4 or higher* is needed with the following
packages:

    - cobra >= 0.18.1
    - requests >= 2.24.0
    - Escher >= 1.7.3

Installation
============

In order to install CobraMod, clone the development branch and install with
using pip with the command ::

    pip install .

By default, required packages will be installed when installing cobramod.

.. note::
    For development, install it with pip and add the argument :code:`-e`. An
    environment YAML file is provided to use it with *conda*.

As mentioned by the original COBRApy developers, while CobraMod is better in a
*conda* environment, it is not necesary. However, an installation in a virtual
environment is encouraged.

How to use it
=============

As with other python packages, to call CobraMod use the import statement in the
current script as :code:`import cobramod`. It is not necesary to import cobra
for CobraMod to work. With CobraMod is is possible to:

- Retrieve data from a database :func:`cobramod.get_data`
- Transform data into COBRApy objects :func:`cobramod.create_object`
- Add reactions from multiple sources :func:`cobramod.add_reactions.`
- Add metabolites from multiple sources :func:`cobramod.add_metabolites`
- Retrieve complete pathways :func:`cobramod.add_pathway`

The posible source for reactions and metabolites are:

- a file with the identifier of the objects, or the whole attributes for a
  custom object.
- a string with the identifier of the object.
- using regular COBRApy objects.

Check the corresponding docstring of the functions for more information or
read the :doc:`guide` to learn about these functions.

To know the database that work in CobraMod, load
:func:`cobramod.available_databases`.
