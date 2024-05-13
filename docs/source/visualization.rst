================
Visualization
================

Cobramod provides two ways to visualization pathways. On the one hand using `Escher <https://github.com/zakandrewking/escher>`_ to
create a 2-dimensional representation and on the other hand using
`3d-force-directed-graph <https://github.com/vasturiano/3d-force-graph?tab=readme-ov-file>`_ to create an interactive three-dimensional representation using threejs and webgl.

----------
 Escher
----------

.. admonition:: Deprecated
   :class: caution

    .. deprecated:: 1.4.0
        The original Python integration of Escher will be removed in a future release due to dependency conflicts with
        Jupyter. It is no longer listed as a dependency in CobraMod and is therefore not automatically installed or
        tested as of version 1.4.0. It is only available to maintain compatibility with older environments.
        The integration embedded in CobraMod will take its place in the future. (See `Escher-custom`_)

CobraMod uses Escher to visualize pathways and fluxes. Each CobraMod pathway includes a visualization method
Pathway.visualize() which automatically generates pathway maps of the respective set of reactions. These
pathway maps can be easily customized to visualize flux distributions using default or user-defined
colors and gradients (linear or quantile normalized).

In the following example, we call the function visualize without any arguments.

.. code-block:: python

    import webbrowser

    pathway.color_negative = "red"
    pathway.color_positive = "green"
    pathway.visualize(solution_fluxes=test_solution)
    webbrowser.open("pathway.html")

We can modify the orientation of our pathway by changing the attribute vertical to True.

.. code-block:: python

    # For flux visualization of the group
    solution =  {
        "GLCpts": -2, "G6PDH2r": -2, "PGL": 0.4, "GND": 1
    }
    # Modifying attributes
    test_model.groups.get_by_id("curated_pathway").visualize(
        solution_fluxes=solution
    )

The visualization method can also be called with the argument solution_fluxes. This argument can be a dictionary with the fluxes of the reactions or a COBRApy Solution. CobraMod assigns colors to the flux values based on the chosen normalizaton method.

By default, the visualization method uses the minimal and maximal values of the solution_flux argument. The user can define the color of the values by changing the attributes color_negative and color_positive. CobraMod creates a linear color gradient with zero flux values colored in grey.

In the following example, we create a dictionary with fluxes and we pass it to the visualization method.

.. code-block:: python

    # Modifying attributes
    test_model.groups.get_by_id("curated_pathway").color_negative = "red"
    test_model.groups.get_by_id("curated_pathway").color_positive = "green"
    test_model.groups.get_by_id("curated_pathway").visualize(
        solution_fluxes=solution
    )

The user can also manually set bounds for color gradients by modifying the CobraMod patway attribute color_min_max. In this example we change the bounds to -10 and 10. This option is useful when comparing different flux simulations.

.. code-block:: python

    # Modifying attributes
    test_model.groups.get_by_id("curated_pathway").color_min_max = [-10, 10]
    test_model.groups.get_by_id("curated_pathway").visualize(
        solution_fluxes=solution
    )

In the next example, we use the default settings for the color_min_max attribute by setting the respective entry to None and change the color gradient to orange and light blue. A list of available colors can be found `here <https://www.w3schools.com/cssref/css_colors.php>`_.

.. code-block:: python

    # New flux with high value
    solution =  {
        "GLCpts": -2, "G6PDH2r": -2, "PGL": 0.4, "GND": 1, "Other": 1000
    }
    # Using defaults
    test_model.groups.get_by_id("curated_pathway").color_min_max = None

    test_model.groups.get_by_id("curated_pathway").color_negative = "orange"
    test_model.groups.get_by_id("curated_pathway").color_positive = "lightskyblue"
    test_model.groups.get_by_id("curated_pathway").visualize(
        solution_fluxes=solution
    )

The user can change the color gradient to a quantile normalization. When choosing this option, the color gradient is determined by the quantiles of thesolution_fluxes argument and not by the minimal und minimal flux values. The user can specify this option by changing the attribute color_quantile to True. This option is useful when the fluxes values vary by several orders of magnitude. For instance, in the previous example, we added a reaction to the dictionary with a flux value of 1000. We can see that the positive colors are pale. Thus, in the next example we change the attribute color_quantile and now the colors are much brighter.

.. code-block:: python

    test_model.groups.get_by_id("curated_pathway").color_quantile = True
    test_model.groups.get_by_id("curated_pathway").visualize(
        solution_fluxes=solution
    )

The user can call the Pathway for a summary of the current attributes.

.. code-block:: python

    test_model.groups.get_by_id("curated_pathway")

CobraMod pathway maps are saved as HTML files with the default name pathway.html. The user can specify the file name with the argument filename. In the following example, we name the file curated_pathway.html.

.. code-block:: python

    test_model.groups.get_by_id("curated_pathway").visualize(
        solution_fluxes=solution, filename = "curated_pathway.html"
    )

.. automethod:: cobramod.core.pathway.Pathway.visualize

testen

------------------
 Escher-custom
------------------
.. versionadded:: 1.4.0

Since the original Escher integration for Jupyter enforces dependencies that conflict with current versions of Jupyter Notebook and Jupyter Lab and can either be installed with significant additional effort or sometimes not at all, CobraMod offers its own integration for Escher.

The frontend is not restricted in its functionality, but the options that can be set from Python are currently limited to those that are also used within CobraMod.

.. jupyter-execute::

    from cobramod import Pathway
    from random import randint
    from cobramod.test import textbook_biocyc_groups

    test_model = textbook_biocyc_groups.copy()

    # Test fluxes
    test_pathway = test_model.groups.get_by_id("test_group")
    test_pathway = Pathway._transform(test_pathway)

    test_solution = {
        reaction.id: randint(-4, 4) for reaction in test_pathway.members
    }
    test_pathway.color_negative = "red"
    test_pathway.color_positive = "green"
    builder = test_pathway.visualize(
        solution_fluxes=test_solution,
        vis="escher-custom",
        never_ask_before_quit= True,
    )
    builder

.. autoclass:: cobramod.visualization.escher.EscherIntegration
    :show-inheritance:
    :members:

.. autoclass:: cobramod.visualization.escher.ReactionScale
    :show-inheritance:
    :members:


----------------------
 Force-directed graph
----------------------
.. versionadded:: 1.4.0

Last but not least, CobraMod offers an integration of 3d-force-directed-graph. This enables a direct representation of a cobrapy group, reaction or cobraMod pathway object in a three-dimensional representation.

.. jupyter-execute::

    from random import randint
    from cobramod.test import textbook_biocyc_groups
    from cobramod.visualization.force_graph import ForceGraphIntegration

    model = textbook_biocyc_groups.copy()
    group = model.groups.get_by_id("test_group")


    solution = {
        reaction.id: randint(-4, 4) for reaction in group.members
    }

    w = ForceGraphIntegration()
    w.model = group
    w.solution = solution
    w

.. autoclass:: cobramod.visualization.force_graph.ForceGraphIntegration

    .. autoproperty:: model
    .. autoproperty:: solution
    .. automethod:: save_layout
    .. automethod:: load_layout


