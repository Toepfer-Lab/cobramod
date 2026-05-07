"""
flux_network.py — Plotly-based 2D metabolic flux network widget for Jupyter

An anywidget that renders an interactive Plotly visualisation of flux
distributions on any COBRA-compatible SBML metabolic model, designed
to be integrated into the CobraMod package.

Usage (Jupyter)
---------------
    from flux_network import FLoV
    import cobra

    model = cobra.io.read_sbml_model("AraCore_v2_1.xml")
    widget = FLoV()
    widget.compartments = "config/aracore.toml"
    widget.model = model
    widget.add_view("FBA", model.optimize())
    widget  # display

The implementation is split across four sibling modules:

* ``_flux_config``   — constants, defaults, ``ViewSpec``, TOML loaders.
* ``_flux_helpers``  — stateless model/flux/colour helpers, CSV loaders.
* ``_flux_layout``   — spring layout, A* routing, hubs, ``HubLine`` /
  ``StationPair``, Plotly trace builder.
* ``_flux_builder``  — ``_BuilderMixin`` with ``_compute_stations``,
  ``_compute_regular_rails``, ``_rebuild_figure`` plus the JSON encoder.
"""

from __future__ import annotations

from contextlib import contextmanager
from importlib import resources
from pathlib import Path
from typing import Optional, Union

import anywidget
import cobra
import networkx as nx
import numpy as np
import pandas as pd
from cobra import Solution
from scipy.spatial import ConvexHull
from traitlets import traitlets

from ._flux_builder import _BuilderMixin
from ._flux_config import ACTIVE_EDGE_EPS, ViewSpec, _parse_config
from ._flux_config import _load_toml as _flux_load_toml
from ._flux_helpers import (
    _emit_progress,
    _is_currency,
    _is_proton,
    _load_ignore_csv,
    _load_single_csv,
    _met_comp,
    _parse_spread,
    _rxn_kind,
)
from ._flux_layout import (
    _fruchterman_reingold,
    _hull_curve_samples,
    _scale_to_region,
)


class FLoV(_BuilderMixin, anywidget.AnyWidget):
    """
    Interactive 2D Plotly-based metabolic flux network widget for Jupyter.

    Renders a bipartite graph (reactions <-> metabolites) with
    compartment-aware layout, flux-coloured edges and nodes,
    view switching, search, and metabolite dragging.

    Usage::

        widget = FLoV()
        widget.add_view("FBA", model.optimize())
        widget.compartments = "config/aracore.toml"  # or dict
        widget.model = cobra.io.read_sbml_model("model.xml")  # single build
        widget  # display in Jupyter
    """

    # anywidget trait synced to frontend
    _figure_json = traitlets.Unicode("{}").tag(sync=True)

    def __init__(
        self,
        drop_no_data: bool = True,
        drop_protons: bool = True,
        drop_biomass: bool = True,
        stoich_flux: bool = False,
        scale_compartments: bool = True,
        hull_tension: float = 0.4,
        density_scale: Optional[float] = None,
        color_modifier: float = 1.0,
        **kwargs,
    ):
        super().__init__(**kwargs)

        # Internal state (must exist before setters run)
        self._graph_built = False
        self._cobra_model: Optional[cobra.Model] = None
        self._cfg: dict = _parse_config({})  # defaults
        self._views: list[ViewSpec] = []

        # Graph state (populated on model set)
        self._G: Optional[nx.Graph] = None
        self._pos: dict[str, tuple[float, float]] = {}
        self._rxn_nodes: list[str] = []
        self._met_nodes: list[str] = []
        self._met_to_rxns: dict[str, list[tuple[str, float]]] = {}
        self._rxn_kind_map: dict[str, str] = {}
        self._rxn_substrates: dict[str, list[str]] = {}
        self._rxn_products: dict[str, list[str]] = {}

        # Configuration — use setters so future assignments behave identically
        # _graph_built is False here so no rebuild is triggered during init
        self.on_msg(self._handle_custom_msg)
        self._drop_no_data = drop_no_data
        self._drop_protons = drop_protons
        self._drop_biomass = drop_biomass
        self._stoich_flux = stoich_flux
        self._scale_compartments = scale_compartments
        self._spread: tuple[float, float, float, float] = (2.0, 2.0, 2.0, 2.0)
        self._radial_spread = 3.0
        self._hull_tension = hull_tension
        self._density_scale = density_scale
        self._color_modifier = color_modifier
        self._ignored_rxns: set[str] = set()
        self._ignored_mets: set[str] = set()
        self._dirty_level: Optional[str] = None
        self._animate_flux = False
        self._anim_min_period = 0.4   # seconds — fastest pulse at peak flux
        self._anim_max_period = 3.0   # seconds — slowest pulse at minimal flux

    # -- Properties --

    @property
    def compartments(self) -> dict[str, dict]:
        """Current compartment configuration."""
        return self._cfg["compartments"]

    @compartments.setter
    def compartments(self, value: Union[str, Path, dict]):
        """Set compartments from TOML file path or pre-parsed dict."""
        if isinstance(value, (str, Path)):
            raw = _flux_load_toml(Path(value))
            self._cfg = _parse_config(raw)
        elif isinstance(value, dict):
            self._cfg = _parse_config(value)
        else:
            raise TypeError(f"Expected str, Path, or dict, got {type(value)}")
        self._mark_dirty("graph")

    @property
    def model(self) -> Optional[cobra.Model]:
        """The COBRA metabolic model."""
        return self._cobra_model

    @model.setter
    def model(self, value: cobra.Model):
        self._cobra_model = value
        self._build_graph()
        if self._views:
            self._rebuild_figure()
        self._dirty_level = None

    @property
    def drop_no_data(self) -> bool:
        """Remove reactions with no flux data in any view."""
        return self._drop_no_data

    @drop_no_data.setter
    def drop_no_data(self, value: bool) -> None:
        self._drop_no_data = value
        self._mark_dirty("figure")

    @property
    def stoich_flux(self) -> bool:
        """Scale edge width by stoichiometric coefficient × flux."""
        return self._stoich_flux

    @stoich_flux.setter
    def stoich_flux(self, value: bool) -> None:
        self._stoich_flux = value
        self._mark_dirty("figure")

    @property
    def hull_tension(self) -> float:
        """Smoothness of compartment hull outlines (0=sharp, 1=very round)."""
        return self._hull_tension

    @hull_tension.setter
    def hull_tension(self, value: float) -> None:
        self._hull_tension = value
        self._mark_dirty("figure")

    @property
    def scale_compartments(self) -> bool:
        """Scale compartment regions proportionally to sqrt(n_reactions)."""
        return self._scale_compartments

    @scale_compartments.setter
    def scale_compartments(self, value: bool) -> None:
        self._scale_compartments = value
        self._mark_dirty("layout")

    @property
    def spread(self) -> Union[float, str]:
        """Scale factor(s) for inter-compartment distances.

        Accepts a plain number (applies uniformly) **or** a string of
        ``<number><direction>`` tokens where direction is one of
        ``u`` (up), ``d`` (down), ``l`` (left), ``r`` (right).
        Any omitted direction defaults to ``1.0`` (no change).

        Examples::

            widget.spread = 1.5          # expand all directions equally
            widget.spread = "2u"         # push compartments apart vertically upward
            widget.spread = "2u1.5r"     # up ×2, right ×1.5, other axes unchanged
            widget.spread = "2u0.8d1.5r" # all three independent
        """
        u, d, l, r = self._spread
        if u == d == l == r:
            return u
        parts = []
        for val, letter in zip((u, d, l, r), "udlr"):
            if val != 1.0:
                parts.append(f"{val:g}{letter}")
        return "".join(parts)

    @spread.setter
    def spread(self, value: Union[float, str]) -> None:
        self._spread = _parse_spread(value)
        self._mark_dirty("layout")

    @property
    def radial_spread(self) -> float:
        """Radial spread factor applied after Fruchterman-Reingold layout.

        Nodes are displaced from their compartment centroid by a factor that
        varies linearly from ``radial_spread`` (at the centroid) down to
        ``1 + 0.2 * (radial_spread - 1)`` (at the outermost node), so central
        nodes are pushed apart the most and edge nodes are barely moved.
        ``1.0`` (default) disables the effect entirely.
        """
        return self._radial_spread

    @radial_spread.setter
    def radial_spread(self, value: float) -> None:
        self._radial_spread = float(value)
        self._mark_dirty("layout")

    @property
    def density_scale(self) -> Optional[float]:
        """Visual density scale in [0, 1] applied to marker sizes, edge widths, and
        colour alpha. ``None`` (default) auto-computes the scale from network size,
        shrinking visuals for models with more than ~80 rendered reactions."""
        return self._density_scale

    @density_scale.setter
    def density_scale(self, value: Optional[float]) -> None:
        self._density_scale = value
        self._mark_dirty("figure")

    @property
    def color_modifier(self) -> float:
        """Scale factor for edge colour saturation. 1.0 = default; >1 more vivid, <1 more faded."""
        return self._color_modifier

    @color_modifier.setter
    def color_modifier(self, value: float) -> None:
        self._color_modifier = value
        self._mark_dirty("figure")

    @contextmanager
    def configure(self):
        """Group configuration assignments without rebuilding.

        Configuration changes only mark the widget dirty.  Trigger one build
        explicitly afterwards by assigning ``model`` or calling ``rebuild()``::

            with widget.configure():
                widget.drop_protons = True
                widget.density_scale = 0.6
            widget.model = model  # single rebuild here
        """
        yield

    def _mark_dirty(self, level: str) -> None:
        """Mark that a future explicit rebuild is needed.

        ``level`` is one of ``'figure'``, ``'layout'``, or ``'graph'``
        (each implies all cheaper levels).  Rebuilds are intentionally not run
        from property setters, so several settings can be changed cheaply before
        assigning ``model`` or calling ``rebuild()`` once.
        """
        if not self._graph_built:
            return
        rank = {"figure": 0, "layout": 1, "graph": 2}
        if level not in rank:
            raise ValueError("level must be one of 'figure', 'layout', or 'graph'")
        if self._dirty_level is None or rank[level] > rank[self._dirty_level]:
            self._dirty_level = level

    def rebuild(self) -> None:
        """Explicitly rebuild the widget after changing configuration or views."""
        if self._cobra_model is None:
            raise ValueError("No model set.")

        level = self._dirty_level
        if level == "graph" or not self._graph_built:
            self._build_graph()
        elif level == "layout":
            self._compute_layout()

        if self._views:
            self._rebuild_figure()
        self._dirty_level = None

    def _resolve_density_scale(self, n_rxns: int) -> float:
        """Return the effective density scale: user value if set, else auto from network size."""
        if self._density_scale is not None:
            return float(np.clip(self._density_scale, 0.0, 1.0))
        return max(0.40, min(1.0, (80.0 / max(80, n_rxns)) ** 0.5))

    def diagnose_reaction(self, reaction_id: str) -> dict:
        """Explain how one reaction and its metabolites are represented.

        This is intended for notebook debugging when a reaction has flux but
        some or all of its metabolites do not appear in the rendered network.
        """
        mdl = self._cobra_model
        if mdl is None:
            raise ValueError("No model set.")
        if reaction_id not in mdl.reactions:
            raise KeyError(f"{reaction_id!r} is not in the model.")
        if not self._graph_built:
            self._build_graph()

        cfg = self._cfg
        rxn = mdl.reactions.get_by_id(reaction_id)
        G = self._G
        flux_by_view = {
            view.label: view.flux_dict.get(reaction_id, np.nan)
            for view in self._views
        }

        met_rxn_count: dict[str, int] = {}
        for rxn_o in mdl.reactions:
            if self._drop_biomass and self._rxn_kind_map.get(rxn_o.id) == "biomass":
                continue
            for met_o in rxn_o.metabolites:
                is_p = _is_proton(met_o, cfg)
                if not _is_currency(met_o, cfg) or (is_p and not self._drop_protons):
                    met_rxn_count[met_o.id] = met_rxn_count.get(met_o.id, 0) + 1

        metabolites = []
        for met, stoich in rxn.metabolites.items():
            is_proton = _is_proton(met, cfg)
            is_currency = _is_currency(met, cfg)
            rendered = bool(G is not None and met.id in G)
            has_edge = bool(G is not None and G.has_edge(reaction_id, met.id))
            reasons = []
            if met.id in self._ignored_mets:
                reasons.append("explicitly ignored")
            if is_currency and not (is_proton and not self._drop_protons):
                reasons.append("filtered as currency")
            if is_proton and self._drop_protons:
                reasons.append("filtered as proton")
            if met_rxn_count.get(met.id, 0) < 2:
                reasons.append("< 2 non-currency/proton connections")
            if not reasons and not rendered:
                reasons.append("not present in rendered graph")
            metabolites.append({
                "id": met.id,
                "name": met.name or met.id,
                "compartment": _met_comp(met, cfg),
                "stoich": float(stoich),
                "connection_count": met_rxn_count.get(met.id, 0),
                "currency": bool(is_currency),
                "proton": bool(is_proton),
                "rendered_metabolite": rendered,
                "rendered_edge": has_edge,
                "hidden_reasons": reasons,
            })

        result = {
            "reaction_id": reaction_id,
            "name": rxn.name or reaction_id,
            "kind": self._rxn_kind_map.get(reaction_id, "unknown"),
            "in_rendered_graph": bool(G is not None and reaction_id in G),
            "flux_by_view": flux_by_view,
            "metabolites": metabolites,
        }
        return result

    @property
    def drop_protons(self) -> bool:
        """Exclude proton (H+) metabolites from the graph."""
        return self._drop_protons

    @drop_protons.setter
    def drop_protons(self, value: bool) -> None:
        self._drop_protons = value
        self._mark_dirty("graph")

    @property
    def drop_biomass(self) -> bool:
        """Exclude biomass reactions from metabolite centroid calculation."""
        return self._drop_biomass

    @drop_biomass.setter
    def drop_biomass(self, value: bool) -> None:
        self._drop_biomass = value
        self._mark_dirty("graph")

    @property
    def ignore(self) -> tuple[set[str], set[str]]:
        """Reaction and metabolite IDs excluded from the graph."""
        return self._ignored_rxns, self._ignored_mets

    @ignore.setter
    def ignore(self, value: Union[str, Path, tuple, None]) -> None:
        if value is None:
            self._ignored_rxns, self._ignored_mets = set(), set()
        elif isinstance(value, (str, Path)):
            self._ignored_rxns, self._ignored_mets = _load_ignore_csv(value)
        elif isinstance(value, tuple) and len(value) == 2:
            self._ignored_rxns, self._ignored_mets = set(value[0]), set(value[1])
        else:
            raise TypeError("ignore expects a file path, a (rxn_set, met_set) tuple, or None")
        self._mark_dirty("graph")

    # -- Public API --

    def add_view(
        self,
        label: str,
        flux_data: Union[Solution, dict, pd.Series, str, Path],
        std_data: Optional[dict] = None,
        pipeline_tag: str = "",
        is_diff: bool = False,
    ) -> None:
        """
        Add a flux view to the visualisation.

        Parameters
        ----------
        label : str
            Display name for this view.
        flux_data : Solution, dict, pd.Series, str, or Path
            Flux values. Accepts cobra.Solution, {rxn_id: flux} dict,
            pd.Series, or path to a CSV file.
        std_data : dict, optional
            Standard deviations {rxn_id: std}.
        pipeline_tag : str, optional
            Tag shown in hover text.
        is_diff : bool, optional
            Mark this view as a computed diff (skipped by the
            auto-diff logic that triggers when two views exist).
        """
        if isinstance(flux_data, Solution):
            fd = {str(k): float(v) for k, v in flux_data.fluxes.items()}
        elif isinstance(flux_data, (pd.Series, dict)):
            fd = {str(k): float(v) for k, v in flux_data.items()}
        elif isinstance(flux_data, (str, Path)):
            s = _load_single_csv(Path(flux_data))
            fd = {str(k): float(v) for k, v in s.items()}
        else:
            raise TypeError(f"Unsupported flux_data type: {type(flux_data)}")

        spec = ViewSpec(
            label=label,
            flux_dict=fd,
            std_dict=std_data or {},
            pipeline_tag=pipeline_tag or label,
            is_diff=is_diff,
        )
        for i, existing in enumerate(self._views):
            if existing.label == label:
                self._views[i] = spec
                break
        else:
            self._views.append(spec)
        if self._graph_built:
            self._mark_dirty("figure")
        if not is_diff:
            self._refresh_auto_diff()

    def clear_views(self) -> None:
        """Remove all flux views."""
        self._views.clear()
        self._figure_json = "{}"
        self._mark_dirty("figure")

    def run_fba(self, label: str = "FBA") -> None:
        """Run FBA on the model and add as a view."""
        if self._cobra_model is None:
            raise ValueError("No model set.")
        _emit_progress("Running FBA...")
        sol = self._cobra_model.optimize()
        if sol.status != "optimal":
            raise RuntimeError(f"FBA status: {sol.status}")
        self.add_view(label, sol)
        self.rebuild()

    def run_pfba(self, label: str = "pFBA") -> None:
        """Run parsimonious FBA on the model and add as a view."""
        if self._cobra_model is None:
            raise ValueError("No model set.")
        _emit_progress("Running pFBA...")
        sol = cobra.flux_analysis.pfba(self._cobra_model)
        self.add_view(label, sol)
        self.rebuild()

    def add_diff_view(
        self,
        label_a: str,
        label_b: str,
        label: str = "Δ",
    ) -> None:
        """Add a view showing the per-reaction flux difference (b − a).

        ``label_a`` and ``label_b`` must already be present as views.
        Reactions only in one of the two contribute their own flux (the
        missing side is treated as 0). Appears as another button in the
        toolbar's Flux view group.

        Called automatically once two non-diff views exist (e.g. after
        ``run_fba()`` then ``run_pfba()``); call directly if you want a
        custom label or a diff between non-adjacent views.
        """
        view_by_label = {v.label: v for v in self._views}
        missing = [l for l in (label_a, label_b) if l not in view_by_label]
        if missing:
            raise ValueError(f"Unknown view label(s): {missing}")
        fa = view_by_label[label_a].flux_dict
        fb = view_by_label[label_b].flux_dict
        diff = {
            r: d
            for r in set(fa) | set(fb)
            if abs(d := fb.get(r, 0.0) - fa.get(r, 0.0)) > ACTIVE_EDGE_EPS
        }
        self.add_view(
            label, diff,
            pipeline_tag=f"{label_b} − {label_a}",
            is_diff=True,
        )

    def _refresh_auto_diff(self) -> None:
        """Maintain a single auto-diff view between the two non-diff views.

        Triggered from ``add_view``: drops any stale auto-diff and rebuilds
        it from the current pair so the diff stays in sync when a view is
        replaced. Skips when there are not exactly two non-diff views.
        """
        non_diff = [v for v in self._views if not v.is_diff]
        self._views = [v for v in self._views if not v.is_diff]
        if len(non_diff) == 2:
            self.add_diff_view(non_diff[0].label, non_diff[1].label)

    def compare_models(
        self,
        model_a: cobra.Model,
        model_b: cobra.Model,
        label_a: str = "A",
        label_b: str = "B",
        method: str = "fba",
    ) -> None:
        """Solve two models and load both as views to diff fluxes.

        Each model is solved independently with its own bounds, then the
        two are merged into a union graph (every reaction from either
        model is drawn) and both solutions are added as side-by-side
        views — analogous to comparing FBA vs pFBA.

        Reactions absent from one model show as no-data in that view, so
        the typical "wildtype vs variant" case overlays cleanly while
        genuinely divergent models still render their full union.

        Parameters
        ----------
        model_a, model_b : cobra.Model
            The two models to compare.
        label_a, label_b : str
            Display names for the two views.
        method : {"fba", "pfba"}
            Flux method to apply to both models.
        """
        if method == "fba":
            solve = lambda m: m.optimize()
        elif method == "pfba":
            solve = lambda m: cobra.flux_analysis.pfba(m)
        else:
            raise ValueError(f"method must be 'fba' or 'pfba', got {method!r}")

        _emit_progress(f"Solving {label_a} ({method})...")
        sol_a = solve(model_a)
        _emit_progress(f"Solving {label_b} ({method})...")
        sol_b = solve(model_b)

        _emit_progress("Merging models for union layout...")
        union = model_a.merge(model_b, inplace=False, objective="left")

        self.model = union
        self.add_view(label_a, sol_a)
        self.add_view(label_b, sol_b)
        self.rebuild()

    def _handle_custom_msg(self, _widget, content: dict, _buffers) -> None:
        if content.get("type") == "run_pfba":
            self.run_pfba()
        elif content.get("type") == "run_fba":
            self.run_fba()

    # -- Internal: graph building --

    def _build_graph(self) -> None:
        """Build the bipartite reaction-metabolite graph from the COBRA model."""
        mdl = self._cobra_model
        if mdl is None:
            return
        _emit_progress("Building graph...")
        cfg = self._cfg

        # Classify reactions
        self._rxn_kind_map = {}
        exch_comp: dict[str, str] = {}
        for rxn in mdl.reactions:
            k = _rxn_kind(rxn, cfg)
            self._rxn_kind_map[rxn.id] = k
            if k != "regular":
                mets = list(rxn.metabolites.keys())
                exch_comp[rxn.id] = _met_comp(mets[0], cfg) if mets else cfg["exchange_key"]

        _kind_counts: dict[str, int] = {}
        for k in self._rxn_kind_map.values():
            _kind_counts[k] = _kind_counts.get(k, 0) + 1
        _emit_progress(
            f"  {len(mdl.reactions)} reactions  "
            f"({_kind_counts.get('regular', 0)} regular, "
            f"{_kind_counts.get('export', 0)} exchange, "
            f"{_kind_counts.get('biomass', 0)} biomass/demand, "
            f"{_kind_counts.get('import', 0)} import)  |  "
            f"{len(mdl.metabolites)} metabolites  |  "
            f"{len(cfg['compartments'])} compartments"
        )

        # Build substrates/products for hover
        self._rxn_substrates = {}
        self._rxn_products = {}
        for rxn in mdl.reactions:
            self._rxn_substrates[rxn.id] = [
                f"{abs(s):.1f}x {m.name or m.id}"
                for m, s in rxn.metabolites.items() if s < 0
            ]
            self._rxn_products[rxn.id] = [
                f"{s:.1f}x {m.name or m.id}"
                for m, s in rxn.metabolites.items() if s > 0
            ]

        # Count metabolite connections (excluding currency).
        # Protons are currency but must be counted when drop_protons=False,
        # otherwise their count stays 0 and the ">= 2" filter discards them.
        met_rxn_count: dict[str, int] = {}
        for rxn in mdl.reactions:
            if self._drop_biomass and self._rxn_kind_map.get(rxn.id) == "biomass":
                continue
            for met in rxn.metabolites:
                is_p = _is_proton(met, cfg)
                if not _is_currency(met, cfg) or (is_p and not self._drop_protons):
                    met_rxn_count[met.id] = met_rxn_count.get(met.id, 0) + 1

        self._met_to_rxns = {}
        for rxn in mdl.reactions:
            if self._drop_biomass and self._rxn_kind_map.get(rxn.id) == "biomass":
                continue
            for met, stoich in rxn.metabolites.items():
                is_proton = _is_proton(met, cfg)
                # When drop_protons=False the user explicitly wants protons
                # visible, so don't let the currency filter swallow them first.
                if _is_currency(met, cfg) and not (is_proton and not self._drop_protons):
                    continue
                if is_proton and self._drop_protons:
                    continue
                if met_rxn_count.get(met.id, 0) < 2:
                    continue
                self._met_to_rxns.setdefault(met.id, []).append(
                    (rxn.id, float(stoich))
                )

        # Build networkx graph
        G = nx.Graph()
        for rxn in mdl.reactions:
            if rxn.id in self._ignored_rxns:
                continue
            if self._drop_biomass and self._rxn_kind_map[rxn.id] == "biomass":
                continue
            if self._rxn_kind_map[rxn.id] == "regular":
                comp_counts: dict[str, int] = {}
                for met_r in rxn.metabolites:
                    mc_r = _met_comp(met_r, cfg)
                    comp_counts[mc_r] = comp_counts.get(mc_r, 0) + 1
                non_exch = {k: v for k, v in comp_counts.items()
                            if k != cfg["exchange_key"]}
                pool = non_exch if non_exch else comp_counts
                comp = max(pool, key=pool.__getitem__) if pool else cfg["exchange_key"]
            elif self._rxn_kind_map[rxn.id] == "transport":
                comp_counts_t: dict[str, int] = {}
                for met_r in rxn.metabolites:
                    mc_r = _met_comp(met_r, cfg)
                    if mc_r != cfg["exchange_key"]:
                        comp_counts_t[mc_r] = comp_counts_t.get(mc_r, 0) + 1
                comp = max(comp_counts_t, key=comp_counts_t.__getitem__) if comp_counts_t else cfg["exchange_key"]
            else:
                comp = exch_comp.get(rxn.id, cfg["exchange_key"])
            G.add_node(rxn.id, ntype="rxn", comp=comp,
                       name=rxn.name or rxn.id, kind=self._rxn_kind_map[rxn.id])
        for met_id, rxn_list in self._met_to_rxns.items():
            if met_id in self._ignored_mets:
                continue
            met_obj = mdl.metabolites.get_by_id(met_id)
            mc = _met_comp(met_obj, cfg)
            G.add_node(met_id, ntype="met", comp=mc, name=met_obj.name or met_id,
                       is_proton=_is_proton(met_obj, cfg))
            for rxn_id, stoich in rxn_list:
                if rxn_id not in self._ignored_rxns and rxn_id in G:
                    G.add_edge(rxn_id, met_id, stoich=stoich)

        # Currency/proton filtering can leave reactions with no rendered
        # metabolite neighbours at all. Drop those orphaned reactions so
        # proton-only or energy-carrier-only fluxes do not appear as lone nodes.
        orphan_rxns = [
            node_id
            for node_id, data in G.nodes(data=True)
            if data["ntype"] == "rxn" and G.degree(node_id) == 0
        ]
        if orphan_rxns:
            G.remove_nodes_from(orphan_rxns)

        self._G = G
        self._rxn_nodes = [n for n, d in G.nodes(data=True) if d["ntype"] == "rxn"]
        self._met_nodes = [n for n, d in G.nodes(data=True) if d["ntype"] == "met"]

        # Layout
        _emit_progress("Computing layout...")
        self._compute_layout()
        self._graph_built = True
        if self._ignored_rxns or self._ignored_mets:
            self._report_ignore(mdl, cfg, met_rxn_count)

    def _report_ignore(
        self, mdl: cobra.Model, cfg: dict, met_rxn_count: dict[str, int]
    ) -> None:
        """Print a summary of ignored IDs, noting if they were already hidden."""
        all_rxn_ids = {r.id for r in mdl.reactions}
        all_met_ids = {m.id for m in mdl.metabolites}

        # Merge flux dicts across all views for zero-flux check
        all_fluxes: dict[str, float] = {}
        for spec in self._views:
            all_fluxes.update(spec.flux_dict)

        if self._ignored_rxns:
            print("── Ignored reactions ──")
            for rid in sorted(self._ignored_rxns):
                if rid not in all_rxn_ids:
                    print(f"  {rid}: not in model (check ID)")
                elif all_fluxes and abs(all_fluxes.get(rid, 0.0)) < 1e-9:
                    print(f"  {rid}: redundant — already carries no flux")
                else:
                    print(f"  {rid}: removed from graph")

        if self._ignored_mets:
            print("── Ignored metabolites ──")
            for mid in sorted(self._ignored_mets):
                if mid not in all_met_ids:
                    print(f"  {mid}: not in model (check ID)")
                    continue
                met_obj = mdl.metabolites.get_by_id(mid)
                if _is_currency(met_obj, cfg):
                    print(f"  {mid}: redundant — already filtered as currency")
                elif _is_proton(met_obj, cfg) and self._drop_protons:
                    print(f"  {mid}: redundant — already filtered as proton")
                elif met_rxn_count.get(mid, 0) < 2:
                    print(f"  {mid}: redundant — already hidden (< 2 connections)")
                else:
                    print(f"  {mid}: removed from graph")

    def _compute_layout(self) -> None:
        """Compute compartment-aware spring layout for all nodes."""
        G = self._G
        if G is None:
            return
        cfg = self._cfg
        compartments = cfg["compartments"]
        pos: dict[str, tuple[float, float]] = {}

        _rxn_counts = {
            ck: sum(1 for r in self._rxn_nodes if G.nodes[r]["comp"] == ck)
            for ck in compartments
        }
        _nonempty_counts = [c for c in _rxn_counts.values() if c > 0]
        _ref_count = float(np.mean(_nonempty_counts)) if _nonempty_counts else 1.0

        # Compute scaled regions, then run a cheap bbox-repulsion pass so that
        # compartments that grow into each other are pushed apart until non-overlapping.
        _active_keys = [ck for ck in compartments if _rxn_counts.get(ck, 0) > 0]

        def _scaled_region(comp_key: str) -> tuple[float, float, float, float]:
            x0, x1, y0, y1 = compartments[comp_key]["region"]
            if not self._scale_compartments:
                return x0, x1, y0, y1
            n = _rxn_counts.get(comp_key, 0)
            scale = float(np.sqrt(max(n, 1) / _ref_count))
            cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
            hw = (x1 - x0) / 2.0 * scale
            hh = (y1 - y0) / 2.0 * scale
            return cx - hw, cx + hw, cy - hh, cy + hh

        # Mutable centers — repulsion moves centers only, half-widths stay fixed
        _centers: dict[str, list[float]] = {}
        _half: dict[str, tuple[float, float]] = {}
        for ck in _active_keys:
            x0, x1, y0, y1 = _scaled_region(ck)
            _centers[ck] = [(x0 + x1) / 2.0, (y0 + y1) / 2.0]
            _half[ck] = ((x1 - x0) / 2.0, (y1 - y0) / 2.0)

        _su, _sd, _sl, _sr = self._spread
        if (_su, _sd, _sl, _sr) != (1.0, 1.0, 1.0, 1.0) and len(_active_keys) > 1:
            gcx = sum(c[0] for c in _centers.values()) / len(_centers)
            gcy = sum(c[1] for c in _centers.values()) / len(_centers)
            for ck in _active_keys:
                dx = _centers[ck][0] - gcx
                dy = _centers[ck][1] - gcy
                # Choose scale per axis based on which side of the centroid
                x_sc = _sr if dx >= 0 else _sl
                y_sc = _su if dy >= 0 else _sd
                _centers[ck][0] = gcx + dx * x_sc
                _centers[ck][1] = gcy + dy * y_sc

        if self._scale_compartments and len(_active_keys) > 1:
            _GAP = 0.3          # minimum clearance between boxes after separation
            for _ in range(50): # converges in far fewer iterations for ~5 compartments
                moved = False
                for i, a in enumerate(_active_keys):
                    for b in _active_keys[i + 1:]:
                        ax, ay = _centers[a]
                        bx, by = _centers[b]
                        ahw, ahh = _half[a]
                        bhw, bhh = _half[b]
                        # Overlap on each axis
                        ox = (ahw + bhw + _GAP) - abs(bx - ax)
                        oy = (ahh + bhh + _GAP) - abs(by - ay)
                        if ox <= 0 or oy <= 0:
                            continue  # already separated
                        # Resolve along the axis of least penetration
                        if ox < oy:
                            shift = ox / 2.0 * (1.0 if bx >= ax else -1.0)
                            _centers[a][0] -= shift
                            _centers[b][0] += shift
                        else:
                            shift = oy / 2.0 * (1.0 if by >= ay else -1.0)
                            _centers[a][1] -= shift
                            _centers[b][1] += shift
                        moved = True
                if not moved:
                    break

        def _effective_region(comp_key: str) -> tuple[float, float, float, float]:
            if comp_key in _centers:
                cx, cy = _centers[comp_key]
                hw, hh = _half[comp_key]
                return cx - hw, cx + hw, cy - hh, cy + hh
            # Inactive compartment — return base region unchanged
            x0, x1, y0, y1 = compartments[comp_key]["region"]
            return x0, x1, y0, y1

        # Per-compartment spring layout for reactions
        for comp_key in compartments:
            x0, x1, y0, y1 = _effective_region(comp_key)
            rxns_here = [r for r in self._rxn_nodes if G.nodes[r]["comp"] == comp_key]
            if not rxns_here:
                continue
            _emit_progress(
                f"Laying out {compartments[comp_key]['label']}..."
            )
            if len(rxns_here) == 1:
                pos[rxns_here[0]] = ((x0 + x1) / 2, (y0 + y1) / 2)
                continue
            rxn_set = set(rxns_here)
            local_G = nx.Graph()
            local_G.add_nodes_from(rxns_here)
            for met_id_l, rxn_list in self._met_to_rxns.items():
                local_rxns = [r for r, _ in rxn_list if r in rxn_set]
                for i in range(len(local_rxns)):
                    for j in range(i + 1, len(local_rxns)):
                        local_G.add_edge(local_rxns[i], local_rxns[j])
            lp = _fruchterman_reingold(
                local_G, k=1.4, iterations=50,
                seed=abs(hash(comp_key)) % (2 ** 31),
            )
            pos.update(_scale_to_region(lp, x0, x1, y0, y1))

        # Metabolite positions — uniform centroid over same-compartment
        # reaction neighbours only. No cross-compartment fallback.
        for i, met_id in enumerate(self._met_nodes):
            met_c = G.nodes[met_id].get("comp", cfg["exchange_key"])
            nbrs = [n for n in G.neighbors(met_id)
                    if G.nodes[n]["ntype"] == "rxn" and n in pos
                    and G.nodes[n].get("comp") == met_c]
            if nbrs:
                pts = np.array([[pos[n][0], pos[n][1]] for n in nbrs], dtype=float)
                pos[met_id] = (float(pts[:, 0].mean()), float(pts[:, 1].mean()))
            elif met_c in compartments:
                x0, x1, y0, y1 = _effective_region(met_c)
                pos[met_id] = ((x0 + x1) / 2.0, (y0 + y1) / 2.0)
            else:
                pos[met_id] = (0.0, 0.0)

        self._pos = pos

        # Radial spread: re-scale each node's distance from its compartment
        # centroid. The scale factor decays linearly from `radial_spread` at
        # the centre to `1 + 0.2*(radial_spread-1)` at the outermost node,
        # so inner nodes are pushed apart more than already-peripheral ones.
        if self._radial_spread != 1.0:
            G = self._G
            for comp_key in compartments:
                nodes_here = [
                    n for n in (self._rxn_nodes + self._met_nodes)
                    if G.nodes[n].get("comp") == comp_key and n in pos
                ]
                if len(nodes_here) < 2:
                    continue
                pts = np.array([[pos[n][0], pos[n][1]] for n in nodes_here], dtype=float)
                cx, cy = pts.mean(axis=0)
                deltas = pts - np.array([cx, cy])
                dists = np.hypot(deltas[:, 0], deltas[:, 1])
                max_dist = dists.max()
                if max_dist == 0.0:
                    continue
                r_norm = dists / max_dist  # 0 = at centroid, 1 = outermost
                f = self._radial_spread
                # scale factor: f at centre, 1+0.2*(f-1) at edge
                scale = 1.0 + (f - 1.0) * (1.0 - 0.99 * r_norm)
                new_pts = np.array([cx, cy]) + deltas * scale[:, np.newaxis]
                for i, n in enumerate(nodes_here):
                    pos[n] = (float(new_pts[i, 0]), float(new_pts[i, 1]))

        # Store effective region function for later use
        self._effective_region = _effective_region

        # ── Compartment hulls (shared by station builder & figure renderer) ──
        # Computed once here so view switches don't re-run the convex-hull pass
        # and stations stay anchored on the same curve as the rendered hulls.
        self._compartment_hulls: dict[str, np.ndarray] = {}
        self._compartment_curve_samples: dict[str, np.ndarray] = {}
        self._compartment_curve_tangents: dict[str, np.ndarray] = {}
        for comp_key, _cfg_c in compartments.items():
            regular_here = [r for r in self._rxn_nodes
                            if G.nodes[r]["comp"] == comp_key
                            and G.nodes[r]["kind"] == "regular"
                            and r in pos]
            if len(regular_here) < 3:
                continue
            pts_h = np.array([[pos[r][0], pos[r][1]] for r in regular_here])
            try:
                hull = ConvexHull(pts_h)
            except Exception:
                continue
            verts = pts_h[hull.vertices]
            centre = verts.mean(axis=0)
            verts = centre + (verts - centre) * 1.22  # match render-time inflation
            self._compartment_hulls[comp_key] = verts
            samples, tangents = _hull_curve_samples(verts, self._hull_tension)
            self._compartment_curve_samples[comp_key] = samples
            self._compartment_curve_tangents[comp_key] = tangents

    # -- Frontend ESM --

    _esm = resources.files("cobramod.static").joinpath("flux_network.mjs").read_text()
