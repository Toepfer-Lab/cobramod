Capabilities
============

With CobraMod is is possible to:

- Retrieve data from a database :func:`cobramod.get_data`
- Transform data into COBRApy objects :func:`cobramod.create_object`
- Add reactions from multiple sources :func:`cobramod.add_reactions`
- Add metabolites from multiple sources :func:`cobramod.add_metabolites`
- Retrieve complete pathways :func:`cobramod.add_pathway`

The posible source for reactions and metabolites are:

- a file with the identifier of the objects, or the whole attributes for a
  custom object.
- a string with the identifier of the object.
- using regular COBRApy objects.

Additionally, CobraMod includes a new type of :class:`cobra.core.group.Group`,
namely :class:`cobramod.Pathway`, which includes new methods such as
:func:`cobramod.Pathway.visualize`

Check the corresponding docstring of the functions for more information or
read the :doc:module/cobramod to learn about these functions.

To know the databases that work in CobraMod, load
:obj:`cobramod.available_databases`.
