"""FLOV — Flux Layout Over Visualisation.

Plotly-based 2D metabolic flux network widget for Jupyter, split across
five modules for maintainability:

* :mod:`._flux_config`   — constants, defaults, ``ViewSpec``, TOML loaders.
* :mod:`._flux_helpers`  — stateless model/flux/colour helpers, CSV loaders.
* :mod:`._flux_layout`   — spring layout, A* routing, ``HubLine`` /
  ``StationPair``, hubs, Plotly trace builder.
* :mod:`._flux_builder`  — ``_BuilderMixin`` with ``_compute_stations``,
  ``_compute_regular_rails``, ``_rebuild_figure`` plus the JSON encoder.
* :mod:`.flux_network`   — the public ``FLoV`` widget class.

The package re-exports :class:`FLoV` so callers can use::

    from cobramod.visualization.FLOV import FLoV
"""

from .flux_network import FLoV

__all__ = ["FLoV"]
