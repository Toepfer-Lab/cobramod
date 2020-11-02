Requirements
============

To run cobramod, *Python 3.7.4 or higher* is needed with the following
packages:

    - cobra >= 0.18.1
    - requests >= 2.24.0

Installation
============

In order to install cobramod, clone the development branch and install with
using pip
with the command ::

    pip install .

By default, required packages will be install when installing cobramod.

.. note::
    For development, install it with pip and add the argument :code:`-e`. An
    environment YAML file is provided to use it with *conda*.

As mentioned by the original COBRApy developers, while cobramod is better in a
*conda* environment, it is not necesary.

How to use it
=============

As with other python packages, to call cobramod use the import statement in the
current script as :code:`import cobramod`. It is not necesary to import cobra
for cobramod to work. With cobramod is is possible to:

- Retrieve data from a database :func:`cobramod.get_data`
- Transform data into COBRApy objects :func:`cobramod.create_object`
- Use files to create reaction :func:`cobramod.add_reactions_from_file`
- Use files to add metabolites :func:`cobramod.add_meta_from_file`
- Retrieve complete pathways :func:`cobramod.add_graph_to_model`
