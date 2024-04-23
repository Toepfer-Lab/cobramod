============
Visualization
============

Cobramod provides two ways to visualization pathways. On the one hand using Escher to
create a 2-dimensional representation and on the other hand using
3d-force-directed-graph to create an interactive three-dimensional representation using threejs and webgl.


In preparation, we create a pathway. In CobraPy, pathways are represented as groups.
However, groups do not have to correspond to pathways, but can be arbitrary combinations of any number of
of any number of metabolites, reactions or other groups.

In the following we will add the pathway with the ID 'SALVADEHYPOX-PWY' from BioCyc using CobraMod if it does not already exist and save it.
from BioCyc and save it as a variable. This pathway will then be created using
Escher and 3d-force-directed-graph.

.. code-block:: python

    from cobramod.test import textbook_biocyc
    from cobramod.visualization.force_graph import ForceGraphBuilder

    test_model = textbook_biocyc.copy()
    add_pathway(
        model=test_model,
        pathway="SALVADEHYPOX-PWY",
        compartment="c",
        directory=dir_data,
        database="ECO",
        ignore_list=[],
        show_imbalance=False,
    )

    pathway = test_model.groups.get_by_id("SALVADEHYPOX-PWY")

    test_solution = {
        reaction.id: randint(-4, 4) for reaction in test_pathway.members
    }

----------
 Escher
----------

.. code-block:: python

    import webbrowser

    pathway.color_negative = "red"
    pathway.color_positive = "green"
    pathway.visualize(solution_fluxes=test_solution)
    webbrowser.open("pathway.html")


.. raw:: html

    <iframe
        src="_static/escher-example.html"
        frameborder="0"
        scrolling="no"
        height= "400"
        width=100%>
    </iframe>

----------
 Escher-custom
----------
Since the original Escher integration for Jupyter enforces dependencies that conflict with current versions of Jupyter Notebook and Jupyter Lab and can either be installed with significant additional effort or sometimes not at all, CobraMod offers its own integration for Escher.

The frontend is not restricted in its functionality, but the options that can be set from Python are currently limited to those that are also used within CobraMod.

----------
 Force-directed graph
----------

Last but not least, CobraMod offers an integration of 3d-force-directed-graph. This enables a direct representation of a cobrapy group, reaction or cobraMod pathway object in a three-dimensional representation.

.. code-block:: python

    from cobramod.visualization.force_graph import ForceGraphBuilder

    builder = ForceGraphBuilder(test_pathway, solution= test_solution)
    builder.open_html()

.. raw:: html

    <iframe
        src="_static/3d-force-example.html"
        frameborder="0"
        scrolling="no"
        height= "400"
        width=100%>
    </iframe>