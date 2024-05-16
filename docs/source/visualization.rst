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
