================
Visualization
================

CobraMod provides users with two different options for visualizing pathways. The first option is
`Escher <https://github.com/zakandrewking/escher>`_, which allows the creation of a 2-dimensional
representation. The second option is
'`3d-force-directed-graph <https://github.com/vasturiano/3d-force-graph?tab=readme-ov-file>`_',
which enables the creation of an interactive three-dimensional representation using
`three.js <https://github.com/mrdoob/three.js/>`_
and
`WebGL <https://www.khronos.org/webgl/>`_.

CobraMod effortlessly combines both tools into the curation and retrieval process. Additionally, they can
operate separately to provide visual representations of `CobraPy <https://opencobra.github.io/cobrapy/>`_ objects.

-----------
 Escher
-----------
.. deprecated:: 1.3.0
    The original Python integration of Escher will be removed in an upcoming release due to dependency
    conflicts with Jupyter. Starting with version 1.3.0 of CobraMod, it will no longer be included as
    a dependency, nor will it be automatically installed or tested. It will only be accessible to ensure
    compatibility with older environments. In the future, the integration incorporated within CobraMod
    will replace it. (Refer to `Escher-custom`_ for more information.)

CobraMod leverages Escher for visualizing pathways and fluxes. Each CobraMod :py:class:`~cobramod.Pathway` comes equipped with the
:py:meth:`~cobramod.Pathway.visualize` method, which effortlessly creates pathway maps for the associated reactions.
These pathway maps can be displayed using `Escher <https://github.com/zakandrewking/escher>`_ and enhanced with flux
distributions. Furthermore, various color schemes and gradients, including linear or quantile normalized options,
can be applied to enrich the visualization.

The following is a simple example that includes flux values and a user-defined color scheme.

.. warning::
    The example may or may not work for you due to dependency issues within Escher's Python binding, as described
    in the deprecation warning.

.. jupyter-input::

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
    )
    builder


.. raw:: html

     <iframe src="./_static/escher-example.html" height="345px" width="100%"></iframe>


.. autoclass:: cobramod.Pathway
    :class-doc-from: class
    :members: vertical, color_negative, color_positive, color_min_max, color_quantile, color_n_steps, color_max_steps, visualize
    :exclude-members: id

------------------
 Escher-custom
------------------
.. versionadded:: 1.3.0

The original `Escher <https://github.com/zakandrewking/escher>`_ integration for `Jupyter <https://jupyter.org/>`_ has dependencies that conflict with current versions of `Jupyter
Notebook <https://jupyter-notebook.readthedocs.io/en/latest/>`_ and `Jupyter Lab <https://jupyterlab.readthedocs.io/en/latest/#>`_. As a result, it can be challenging to install or may not work at all. Consequently,
installation can be pretty challenging or may not work altogether. CobraMod offers its own Escher integration,
which currently has limitations on the options that can be set from Python, although the frontend functionality
remains unrestricted.

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
.. versionadded:: 1.3.0

CobraMod provides an integration of
'`3d-force-directed-graph <https://github.com/vasturiano/3d-force-graph?tab=readme-ov-file>`_'
that enables a direct representation of a
CobraPy :py:class:`~cobra.core.group.Group`, :py:class:`~cobra.Reaction`, or CobraMod :py:class:`~cobramod.Pathway` object in a three-dimensional format. This offers the
benefit of increased space compared to a two-dimensional representation. Moreover, it facilitates the
visualization of data, such as fluxes, through animated representations, thereby enhancing the overall
comprehension of the information.

The same pathway previously visualized using Escher is now presented in three dimensions below. Additional
customization options can be found in the class definition of the :py:class:`~cobramod.visualization.force_graph.ForceGraphIntegration`.

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

----------------------
 FloV Graph
----------------------
FLoV is a flux-centred network visualization that draws metabolic models in the style of a railway map.
Reactions and metabolites are arranged into compartment regions, while transport reactions are collected into
station-like hubs and routed as bundled lines between compartments.

The layout starts from a force-directed graph, maps the result onto a compartment grid, and then uses A*
pathfinding to connect metabolites and reactions across compartment boundaries. Most of the graph construction
is handled automatically in the backend; the user mainly provides flux data and, when needed, a small TOML
annotation file describing compartment labels, colours, and regions.

.. list-table::
    :widths: 1 1
    :header-rows: 1

    * - Notebook usage
      - Compartment annotation
    * - .. code-block:: python

            import cobra
            from cobramod.visualization.flux_network import FLoV

            model = cobra.io.read_sbml_model("model.xml")

            widget = FLoV()
            widget.compartments = "config/model_compartments.toml"
            widget.add_view("FBA", model.optimize())
            widget.model = model
            widget
      - .. code-block:: toml

            [compartment.d]
            label = "Stroma"
            region = [-11.5, -3.0, -4.5, 6.0]
            color = "#1e8449"
            fill = "rgba(30,132,73,0.06)"

The ``region`` entry uses ``[x_min, x_max, y_min, y_max]`` coordinates. Adjust these values to move or resize a
compartment in the final graph. The TOML file can also define import, export, and biomass prefixes, plus currency
metabolite keywords, when the model naming conventions differ from the defaults.

Comparing two SBML models
^^^^^^^^^^^^^^^^^^^^^^^^^

FLoV can also compare two SBML models directly. This is useful for checking how a curated model, mutant model,
or alternative reconstruction changes the flux distribution relative to a reference model.

The comparison workflow has three steps:

* each SBML file is loaded as a separate :class:`cobra.Model`;
* both models are solved independently with the selected flux method, either ``"fba"`` or ``"pfba"``;
* FLoV asks COBRApy to merge the two models into a union model, so reactions from either model can be drawn in
  one shared layout.

The two solutions are then added as separate views. A third difference view is created automatically and shows
``label_b - label_a`` for reactions that changed. Reactions that exist only in one model remain visible in the
union graph; the missing side is treated as no data for that view.

.. code-block:: python

    import cobra
    from cobramod.visualization.flux_network import FLoV

    reference = cobra.io.read_sbml_model("models/reference.xml")
    variant = cobra.io.read_sbml_model("models/variant.xml")

    widget = FLoV()
    widget.compartments = "config/model_compartments.toml"
    widget.compare_models(
        reference,
        variant,
        label_a="Reference",
        label_b="Variant",
        method="pfba",
    )
    widget

Internally, the merged layout is built with COBRApy's model merge operation. The first model supplies the left
objective during the merge, while the flux values shown in the views still come from the two independently solved
input models. This keeps the visual comparison focused on flux changes instead of forcing both models through a
single shared optimization problem.

Configuration flags
^^^^^^^^^^^^^^^^^^^

FLoV can be configured either when the widget is created or by assigning properties before the model is loaded.
The spacing flags are especially important because they determine whether compartment regions and station routes
have enough room to separate cleanly. Settings that affect graph membership, layout, or styling are marked dirty
instead of rebuilding immediately, so several options can be changed together and applied once.

.. code-block:: python

    from cobra.flux_analysis import pfba
    from cobramod.visualization.flux_network import FLoV

    TOML_PATH = "TOMLPATH"
    fba_sol = model.optimize()
    pfba_sol = pfba(model)

    widget = FLoV()
    with widget.configure():
        widget.compartments = TOML_PATH  
        widget.density_scale = 0.7
        widget.spread = 5
        widget.radial_spread = 8.0
        widget.drop_protons = True
        widget.scale_compartments = True
        widget.drop_no_data = True
        widget.report_ignore = False
        widget.add_view("pFBA", pfba_sol)
        widget.add_view("FBA", fba_sol)

    widget.model = model
    widget

.. list-table::
    :widths: 1 1 3
    :header-rows: 1

    * - Flag
      - Default
      - Effect
    * - ``drop_no_data``
      - ``True``
      - Removes reactions that have no flux data in any loaded view. Turn this off when you want to inspect the
        full model topology even if many reactions are inactive or missing from the flux input.
    * - ``drop_protons``
      - ``True``
      - Excludes proton metabolites from the rendered graph. Set to ``False`` when proton transfer itself is part
        of the question being inspected.
    * - ``drop_biomass``
      - ``True``
      - Excludes biomass and demand-like reactions from graph construction. This keeps central biomass reactions
        from dominating the layout.
    * - ``stoich_flux``
      - ``False``
      - Scales route widths by ``abs(stoichiometric coefficient * flux)`` instead of flux alone. This highlights
        reactions where large stoichiometric coefficients make a route visually important.
    * - ``scale_compartments``
      - ``True``
      - Expands compartment regions according to the number of rendered reactions assigned to each compartment.
        Disable it when the TOML ``region`` boxes should be used exactly as written.
    * - ``spread``
      - ``2.0``
      - Controls spacing between compartment regions. Use a number for uniform spacing, or directional tokens
        such as ``"2u"``, ``"1.5r"``, or ``"2u0.8d1.5r"``. Directions are ``u`` for up, ``d`` for down,
        ``l`` for left, and ``r`` for right.
    * - ``radial_spread``
      - ``3.0``
      - Pushes nodes away from the centre of their compartment after the force-directed layout. ``1.0`` disables
        the extra radial push; larger values spread crowded centres.
    * - ``hull_tension``
      - ``0.4``
      - Controls how smooth compartment hull outlines are. Lower values keep sharper outlines; higher values make
        hulls rounder. This is mainly a visual tuning option and usually does not need to be changed.
    * - ``density_scale``
      - ``None``
      - Controls marker sizes, edge widths, and colour opacity in the range ``0`` to ``1``. ``None`` lets FLoV
        shrink dense networks automatically. Change this only when the default visual density is too busy or too
        sparse for a specific figure.
    * - ``report_ignore``
      - ``False``
      - Prints a summary of explicitly ignored reaction and metabolite IDs during rebuilds. Enable this when
        debugging an ignore list or checking whether an ID was already hidden by another filter.

Additional inputs can be changed before loading the model, or followed by ``widget.rebuild()`` when a model is
already assigned:

.. list-table::
    :widths: 1 3
    :header-rows: 1

    * - Property
      - Accepted values
    * - ``compartments``
      - TOML file path, :class:`pathlib.Path`, or a parsed dictionary. This sets compartment labels, colours,
        hull fills, regions, exchange/import/export/biomass prefixes, and currency metabolite keywords.
    * - ``ignore``
      - ``None``, a CSV file path, or a ``(reaction_ids, metabolite_ids)`` tuple. Ignored IDs are removed from
        graph construction.
    * - ``diagnose_reaction("RXN_ID")``
      - Debug helper for one reaction ID. It returns a dictionary with the reaction kind, flux values per view,
        whether the reaction appears in the rendered graph, and why individual metabolites may be hidden
        (currency filter, proton filter, ignore list, or too few graph connections). This works seperately no rebuilding.

When changing these properties after a model has already been assigned, call ``widget.rebuild()`` to apply them.
For several edits at once, use ``with widget.configure():`` to group assignments and then rebuild once.

.. autoclass:: cobramod.visualization.flux_network.FLoV
