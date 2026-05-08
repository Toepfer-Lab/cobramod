"""Figure-building functions for the FLoV widget.

Turns a built graph + views into the JSON figure descriptor consumed by the
frontend. Three top-level entry points, each split into one helper per
numbered comment block:

* :func:`compute_stations` — bucket non-regular reactions into station pairs
  and route their spines through the exchange region with A*.
* :func:`compute_regular_rails` — route in-compartment reaction↔metabolite
  edges on a per-compartment grid, assigning lane slots.
* :func:`rebuild_figure` — assemble per-view trace data and serialise to JSON.

State that used to live on ``self`` is bundled into :class:`BuilderContext`,
which the host widget builds before calling these functions.
"""
from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Optional

import cobra
import networkx as nx
import numpy as np

from ._flux_config import (
    ACTIVE_EDGE_EPS,
    MAX_LINES_PER_PAIR,
    MAX_RXNS_PER_LINE,
    ViewSpec,
    _LINE_KLASS_COLORS,
    _LINE_KLASS_INACTIVE,
    _LINE_KLASSES,
    _RAIL_ROUTE_PALETTE,
)
from ._flux_helpers import (
    _edge_bucket,
    _flux_to_hex,
    _has_flux_data,
    _met_class,
    _met_comp,
    _modulate_route_color,
    _pair_key,
)
from ._flux_layout import (
    HubLine,
    StationPair,
    _a_star_grid,
    _grid_snap,
    _nearest_open_cell,
    _polygon_mask,
    _project_to_curve,
    _simplify_polyline,
)


# ==========================================================================
# JSON SERIALISATION
# ==========================================================================

# Convert NumPy values into JSON-serialisable Python objects.
def _json_default(obj):
    """``json.dumps(default=...)`` callable handling numpy scalars/arrays."""
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serialisable")


# ==========================================================================
# CONTEXT + GRID DATACLASSES
# ==========================================================================

@dataclass(frozen=True)
class BuilderContext:
    """Read-only snapshot of FLoV state needed to build the figure JSON."""
    cobra_model: Optional[cobra.Model]
    cfg: dict
    G: Optional[nx.Graph]
    pos: dict[str, tuple[float, float]]
    views: list[ViewSpec]
    rxn_nodes: list[str]
    met_nodes: list[str]
    met_to_rxns: dict[str, list[tuple[str, float]]]
    rxn_kind_map: dict[str, str]
    rxn_substrates: dict[str, list[str]]
    rxn_products: dict[str, list[str]]
    compartment_hulls: dict[str, np.ndarray]
    compartment_curve_samples: dict[str, np.ndarray]
    compartment_curve_tangents: dict[str, np.ndarray]
    drop_no_data: bool
    stoich_flux: bool
    hull_tension: float
    effective_region: Callable[[str], tuple[float, float, float, float]]
    resolve_density_scale: Callable[[int], float]


@dataclass(frozen=True)
class RouteGrid:
    """Forbidden-cell routing grid covering the global station-routing area."""
    x_min: float
    y_min: float
    step: float
    n_rows: int
    n_cols: int

    # Convert an x/y coordinate into a clamped grid cell.
    def to_cell(self, p: tuple[float, float]) -> tuple[int, int]:
        r = int(round((p[1] - self.y_min) / self.step))
        c = int(round((p[0] - self.x_min) / self.step))
        r = max(0, min(self.n_rows - 1, r))
        c = max(0, min(self.n_cols - 1, c))
        return r, c

    # Convert a grid cell back into an x/y coordinate.
    def to_xy(self, cell: tuple[int, int]) -> tuple[float, float]:
        return (self.x_min + cell[1] * self.step,
                self.y_min + cell[0] * self.step)


@dataclass
class CompartmentGrid:
    """Per-compartment routing grid for regular rxn↔met edges."""
    comp: str
    x_min: float
    y_min: float
    step: float
    n_rows: int
    n_cols: int
    allowed: np.ndarray
    allowed_cells: np.ndarray
    allowed_xy: np.ndarray
    blocked: np.ndarray

    # Convert a compartment grid cell back into an x/y coordinate.
    def to_xy(self, cell: tuple[int, int]) -> tuple[float, float]:
        return (
            float(self.x_min + cell[1] * self.step),
            float(self.y_min + cell[0] * self.step),
        )


# ==========================================================================
# STATION HELPERS
# ==========================================================================

# Group active non-regular reactions by station endpoint pair and line class.
def _bucket_pair_reactions(
    ctx: BuilderContext, active_rxn_set: set[str],
) -> dict[tuple[str, str], dict[str, list[str]]]:
    """Section 1: bucket non-regular reactions by (compA, compB) → klass → ids."""
    mdl = ctx.cobra_model
    cfg = ctx.cfg
    exch = cfg["exchange_key"]
    pair_buckets: dict[tuple[str, str], dict[str, list[str]]] = {}
    for rxn in mdl.reactions:
        if rxn.id not in active_rxn_set:
            continue
        kind = ctx.rxn_kind_map.get(rxn.id, "regular")
        if kind not in ("transport", "import", "export", "biomass"):
            continue
        non_exch = {_met_comp(m, cfg) for m in rxn.metabolites} - {exch}
        pk = _pair_key(kind, non_exch, cfg)
        if pk is None:
            continue
        if pk[0] not in ctx.compartment_hulls and pk[0] != exch:
            continue
        if pk[1] not in ctx.compartment_hulls and pk[1] != exch:
            continue
        klass = _met_class(rxn, cfg)
        pair_buckets.setdefault(pk, {}).setdefault(klass, []).append(rxn.id)
    return pair_buckets


# Turn pair buckets into ordered station line objects.
def _build_station_lines(
    pair_buckets: dict[tuple[str, str], dict[str, list[str]]],
) -> list[StationPair]:
    """Section 2: assemble StationPair objects with ordered HubLines."""
    stations: list[StationPair] = []
    for (a, b), klass_map in pair_buckets.items():
        ordered_lines: list[HubLine] = []
        for ki, klass in enumerate(_LINE_KLASSES):
            rxn_ids = klass_map.get(klass, [])
            if not rxn_ids:
                continue
            if MAX_RXNS_PER_LINE is not None:
                rxn_ids = rxn_ids[:MAX_RXNS_PER_LINE]
            ordered_lines.append(HubLine(
                klass=klass,
                rxn_ids=list(rxn_ids),
                spine_offset_idx=ki,
            ))
            if len(ordered_lines) >= MAX_LINES_PER_PAIR:
                break
        if not ordered_lines:
            continue
        stations.append(StationPair(
            pair_id=f"{a}|{b}",
            comp_a=a, comp_b=b,
            lines=ordered_lines,
        ))
    return stations


# Compute the padded global routing bounds and grid spacing.
def _compute_routing_grid_extent(
    ctx: BuilderContext,
) -> tuple[float, float, float, float, float]:
    """Section 3 prelude: padded x/y bounds + grid step from current layout."""
    pos = ctx.pos
    all_x = [p[0] for p in pos.values()]
    all_y = [p[1] for p in pos.values()]
    x_min, x_max = (min(all_x), max(all_x)) if all_x else (-12.0, 8.0)
    y_min, y_max = (min(all_y), max(all_y)) if all_y else (-9.0, 7.0)
    pad = 2.0
    x_min -= pad
    x_max += pad
    y_min -= pad
    y_max += pad
    span = max(x_max - x_min, y_max - y_min)
    step = max(0.4, span / 50.0)
    return x_min, x_max, y_min, y_max, step


# Resolve the centroid used to aim station anchors for a compartment.
def _comp_centroid(
    ctx: BuilderContext,
    ck: str,
    extent: tuple[float, float, float, float, float],
) -> tuple[float, float]:
    x_min, x_max, y_min, y_max, _ = extent
    if ck in ctx.compartment_hulls:
        v = ctx.compartment_hulls[ck]
        return float(v[:, 0].mean()), float(v[:, 1].mean())
    compartments = ctx.cfg["compartments"]
    if ck in compartments:
        x0, x1, y0, y1 = compartments[ck]["region"]
        return (x0 + x1) / 2.0, (y0 + y1) / 2.0
    return (x_min + x_max) / 2.0, (y_min + y_max) / 2.0


# Place a synthetic exchange anchor on the nearest global boundary.
def _exch_anchor(
    other_centroid: tuple[float, float],
    extent: tuple[float, float, float, float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Synthesise an exchange anchor along the bounding-box edge nearest *other_centroid*."""
    x_min, x_max, y_min, y_max, _ = extent
    ox, oy = other_centroid
    cx, cy = (x_min + x_max) / 2.0, (y_min + y_max) / 2.0
    dx, dy = ox - cx, oy - cy
    if abs(dx) > abs(dy):
        ax = x_max - 0.5 if dx > 0 else x_min + 0.5
        ay = oy
        tan = (0.0, 1.0)
    else:
        ay = y_max - 0.5 if dy > 0 else y_min + 0.5
        ax = ox
        tan = (1.0, 0.0)
    return (ax, ay), tan


# Project station endpoints onto compartment hull curves.
def _anchor_stations_to_hulls(
    ctx: BuilderContext,
    stations: list[StationPair],
    extent: tuple[float, float, float, float, float],
) -> None:
    """Section 3: anchor each station on its compartment's smoothed hull curve."""
    step = extent[4]
    for st in stations:
        cb_x, cb_y = _comp_centroid(ctx, st.comp_b, extent)
        ca_x, ca_y = _comp_centroid(ctx, st.comp_a, extent)
        # Anchor A
        if st.comp_a in ctx.compartment_curve_samples:
            samples = ctx.compartment_curve_samples[st.comp_a]
            tans = ctx.compartment_curve_tangents[st.comp_a]
            idx = _project_to_curve((cb_x, cb_y), samples)
            ax, ay = samples[idx]
            tx, ty = tans[idx]
        else:
            (ax, ay), (tx, ty) = _exch_anchor((ca_x, ca_y), extent)
        sx, sy = _grid_snap(ax, ay, step)
        st.anchor_a = (float(sx), float(sy))
        st.tangent_a = (float(tx), float(ty))
        # Anchor B
        if st.comp_b in ctx.compartment_curve_samples:
            samples = ctx.compartment_curve_samples[st.comp_b]
            tans = ctx.compartment_curve_tangents[st.comp_b]
            idx = _project_to_curve((ca_x, ca_y), samples)
            bx, by = samples[idx]
            tx2, ty2 = tans[idx]
        else:
            (bx, by), (tx2, ty2) = _exch_anchor((cb_x, cb_y), extent)
        sx, sy = _grid_snap(bx, by, step)
        st.anchor_b = (float(sx), float(sy))
        st.tangent_b = (float(tx2), float(ty2))


# Spread nearby station anchors along each compartment curve.
def _deoverlap_station_anchors(
    ctx: BuilderContext, stations: list[StationPair], step: float,
) -> None:
    """Section 4: nudge co-located anchors apart along each compartment curve."""
    for ck, samples in ctx.compartment_curve_samples.items():
        stationed = [(st, "a") for st in stations if st.comp_a == ck] + \
                    [(st, "b") for st in stations if st.comp_b == ck]
        if len(stationed) < 2:
            continue
        indexed = []
        for st, side in stationed:
            anchor = st.anchor_a if side == "a" else st.anchor_b
            idx = _project_to_curve(anchor, samples)
            indexed.append((idx, st, side))
        indexed.sort(key=lambda t: t[0])
        min_gap = max(2, len(samples) // 30)
        prev_idx = -10 ** 9
        for k, (idx, st, side) in enumerate(indexed):
            if idx - prev_idx < min_gap:
                idx = prev_idx + min_gap
                if idx >= len(samples):
                    idx = len(samples) - 1
                pt = samples[idx]
                tan = ctx.compartment_curve_tangents[ck][idx]
                sx, sy = _grid_snap(pt[0], pt[1], step)
                if side == "a":
                    st.anchor_a = (float(sx), float(sy))
                    st.tangent_a = (float(tan[0]), float(tan[1]))
                else:
                    st.anchor_b = (float(sx), float(sy))
                    st.tangent_b = (float(tan[0]), float(tan[1]))
            prev_idx = idx


# Build the forbidden-cell bitmap used by station routing.
def _build_forbidden_grid(
    ctx: BuilderContext,
    active_compartments: set[str],
    extent: tuple[float, float, float, float, float],
) -> tuple[np.ndarray, RouteGrid]:
    """Section 5: build forbidden bitmap from compartment hulls."""
    x_min, x_max, y_min, y_max, step = extent
    n_cols = int(np.ceil((x_max - x_min) / step)) + 1
    n_rows = int(np.ceil((y_max - y_min) / step)) + 1
    cols = np.arange(n_cols)
    rows = np.arange(n_rows)
    gx, gy = np.meshgrid(x_min + cols * step, y_min + rows * step, indexing="xy")
    forbidden = np.zeros((n_rows, n_cols), dtype=bool)
    for ck, verts in ctx.compartment_hulls.items():
        if ctx.drop_no_data and ck not in active_compartments:
            continue
        forbidden |= _polygon_mask(verts, gx, gy)
    grid = RouteGrid(
        x_min=x_min, y_min=y_min, step=step, n_rows=n_rows, n_cols=n_cols,
    )
    return forbidden, grid


# Generate nearby hull anchor candidates for route optimisation.
def _anchor_candidates(
    ctx: BuilderContext,
    comp_key: str,
    current_anchor: tuple[float, float],
    current_tangent: tuple[float, float],
    toward: tuple[float, float],
    step: float,
    max_candidates: int = 14,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    if comp_key not in ctx.compartment_curve_samples:
        return [(current_anchor, current_tangent)]
    samples = ctx.compartment_curve_samples[comp_key]
    tangents = ctx.compartment_curve_tangents[comp_key]
    if len(samples) == 0:
        return [(current_anchor, current_tangent)]

    target_idx = _project_to_curve(toward, samples)
    if target_idx < 0:
        return [(current_anchor, current_tangent)]
    stride = max(1, len(samples) // max_candidates)
    offsets = [0]
    for k in range(1, max_candidates // 2 + 1):
        offsets.extend([k * stride, -k * stride])

    candidates: list[tuple[tuple[float, float], tuple[float, float]]] = []
    seen: set[tuple[float, float]] = set()
    for off in offsets:
        idx = (target_idx + off) % len(samples)
        pt = samples[idx]
        tan = tangents[idx]
        anchor = _grid_snap(float(pt[0]), float(pt[1]), step)
        key = (float(anchor[0]), float(anchor[1]))
        if key in seen:
            continue
        seen.add(key)
        candidates.append((
            (float(anchor[0]), float(anchor[1])),
            (float(tan[0]), float(tan[1])),
        ))
    return candidates or [(current_anchor, current_tangent)]


# Route between two coordinates through the station grid.
def _route_cells_between(
    grid: RouteGrid,
    forbidden: np.ndarray,
    a: tuple[float, float],
    b: tuple[float, float],
) -> Optional[list[tuple[int, int]]]:
    ca = grid.to_cell(a)
    cb = grid.to_cell(b)
    saved = bool(forbidden[ca]), bool(forbidden[cb])
    forbidden[ca] = False
    forbidden[cb] = False
    cells = _a_star_grid(ca, cb, forbidden)
    forbidden[ca] = saved[0]
    forbidden[cb] = saved[1]
    return cells


# Score a routed cell path by total movement length.
def _route_score(cells: Optional[list[tuple[int, int]]]) -> float:
    if not cells:
        return float("inf")
    score = 0.0
    for a, b in zip(cells, cells[1:]):
        dr = a[0] - b[0]
        dc = a[1] - b[1]
        score += float((dr * dr + dc * dc) ** 0.5)
    return score


# Pick the best pair of station anchors for a routed spine.
def _optimise_station_anchors(
    ctx: BuilderContext,
    st: StationPair,
    grid: RouteGrid,
    forbidden: np.ndarray,
    extent: tuple[float, float, float, float, float],
) -> Optional[list[tuple[int, int]]]:
    ca = _comp_centroid(ctx, st.comp_a, extent)
    cb = _comp_centroid(ctx, st.comp_b, extent)
    cand_a = _anchor_candidates(ctx, st.comp_a, st.anchor_a, st.tangent_a, cb, grid.step)
    cand_b = _anchor_candidates(ctx, st.comp_b, st.anchor_b, st.tangent_b, ca, grid.step)
    best: Optional[tuple[
        float, tuple[float, float], tuple[float, float],
        tuple[float, float], tuple[float, float],
        Optional[list[tuple[int, int]]],
    ]] = None
    for anchor_a, tangent_a in cand_a:
        for anchor_b, tangent_b in cand_b:
            cells = _route_cells_between(grid, forbidden, anchor_a, anchor_b)
            score = _route_score(cells)
            if best is None or score < best[0]:
                best = (score, anchor_a, tangent_a, anchor_b, tangent_b, cells)
    if best is None or best[5] is None:
        return None
    st.anchor_a = best[1]
    st.tangent_a = best[2]
    st.anchor_b = best[3]
    st.tangent_b = best[4]
    return best[5]


# Route station spines through the exchange region.
def _route_station_spines(
    ctx: BuilderContext,
    stations: list[StationPair],
    grid: RouteGrid,
    forbidden: np.ndarray,
    extent: tuple[float, float, float, float, float],
) -> None:
    """Section 6: route every station spine via A*, with anchor optimisation."""
    for st in stations:
        cells = _optimise_station_anchors(ctx, st, grid, forbidden, extent)
        if cells is None:
            cells = _route_cells_between(grid, forbidden, st.anchor_a, st.anchor_b)
        if not cells:
            st.spine = [st.anchor_a, st.anchor_b]
            continue
        poly = [grid.to_xy(c) for c in cells]
        poly[0] = st.anchor_a
        poly[-1] = st.anchor_b
        st.spine = _simplify_polyline(poly)


# Collect hidden non-regular reaction-to-metabolite edges for toggling.
def _collect_inner_edges(
    ctx: BuilderContext, active_rxn_set: set[str],
) -> list[dict]:
    """Section 7: rxn→met edges of non-regular reactions (kept for frontend toggling)."""
    G = ctx.G
    pos = ctx.pos
    inner_edges: list[dict] = []
    non_regular_active = {
        r for r in active_rxn_set
        if ctx.rxn_kind_map.get(r) in ("transport", "import", "export", "biomass")
        and r in pos
    }
    for rxn_id in non_regular_active:
        rx, ry = pos[rxn_id]
        for n1, n2 in G.edges(rxn_id):
            met_id = n2 if n1 == rxn_id else n1
            if G.nodes[met_id].get("ntype") != "met":
                continue
            if met_id not in pos:
                continue
            mx, my = pos[met_id]
            inner_edges.append({
                "rxn_id": rxn_id,
                "met_id": met_id,
                "x": [float(rx), float(mx)],
                "y": [float(ry), float(my)],
            })
    return inner_edges


# Build all station geometry for active non-regular reactions.
def compute_stations(
    ctx: BuilderContext,
    active_rxn_set: set[str],
    active_compartments: set[str],
) -> tuple[list[StationPair], float, list[dict]]:
    """Build station list and run A* routing.

    Returns ``(stations, grid_step, inner_edges)``.
    """
    if ctx.cobra_model is None or ctx.G is None or not ctx.pos:
        return [], 0.5, []

    pair_buckets = _bucket_pair_reactions(ctx, active_rxn_set)
    if not pair_buckets:
        return [], 0.5, []

    stations = _build_station_lines(pair_buckets)

    extent = _compute_routing_grid_extent(ctx)
    step = extent[4]

    _anchor_stations_to_hulls(ctx, stations, extent)
    _deoverlap_station_anchors(ctx, stations, step)
    forbidden, grid = _build_forbidden_grid(ctx, active_compartments, extent)
    _route_station_spines(ctx, stations, grid, forbidden, extent)
    inner_edges = _collect_inner_edges(ctx, active_rxn_set)

    return stations, step, inner_edges


# ==========================================================================
# RAIL HELPERS
# ==========================================================================

# Return fallback rail data when no regular routes can be built.
def _default_rail_data(
    ctx: BuilderContext, rxn_nodes: list[str], met_nodes: list[str],
) -> dict:
    """Empty-graph fallback used when there are no regular edges to route."""
    pos = ctx.pos
    return {
        "routes": [],
        "guides": [],
        "rxn_pos": {r: pos.get(r, (0.0, 0.0)) for r in rxn_nodes},
        "met_pos": {m: pos.get(m, (0.0, 0.0)) for m in met_nodes},
        "routed_rxns": set(),
        "routed_mets": set(),
        "grid_step": 0.4,
    }


# Rank reactions by total absolute flux across views.
def _compute_flux_priority(views: list[ViewSpec]) -> dict[str, float]:
    """Sum of |flux| across views per reaction — used to order route assignment."""
    flux_priority: dict[str, float] = defaultdict(float)
    for view in views:
        for rxn_id, flux in view.flux_dict.items():
            if _has_flux_data(float(flux)):
                flux_priority[rxn_id] += abs(float(flux))
    return flux_priority


# Gather regular reaction-metabolite edge records per compartment.
def _collect_edge_records(
    ctx: BuilderContext,
    rxn_nodes: list[str],
    met_nodes: list[str],
    flux_priority: dict[str, float],
) -> tuple[dict[str, list[dict]], dict[str, int]]:
    """Collect per-compartment records of regular rxn↔met edges."""
    G = ctx.G
    edge_records: dict[str, list[dict]] = defaultdict(list)
    node_degree: dict[str, int] = defaultdict(int)
    rxn_set = set(rxn_nodes)
    met_set = set(met_nodes)
    for n1, n2, edata in G.edges(data=True):
        if G.nodes[n1]["ntype"] == "rxn" and G.nodes[n2]["ntype"] == "met":
            rxn_id, met_id = n1, n2
        elif G.nodes[n2]["ntype"] == "rxn" and G.nodes[n1]["ntype"] == "met":
            rxn_id, met_id = n2, n1
        else:
            continue
        if rxn_id not in rxn_set or met_id not in met_set:
            continue
        if ctx.rxn_kind_map.get(rxn_id, "regular") != "regular":
            continue
        comp = G.nodes[rxn_id].get("comp")
        if comp != G.nodes[met_id].get("comp") or comp is None:
            continue
        stoich = float(edata.get("stoich", 1.0))
        line_id = f"{rxn_id}->{met_id}"
        edge_records[comp].append({
            "line_id": line_id,
            "rxn_id": rxn_id,
            "met_id": met_id,
            "comp": comp,
            "stoich": stoich,
            "priority": float(flux_priority.get(rxn_id, 0.0)) + abs(stoich),
        })
        node_degree[rxn_id] += 1
        node_degree[met_id] += 1
    return edge_records, node_degree


# Build the routing grid for one compartment.
def _build_compartment_grid(
    ctx: BuilderContext, comp: str, n_node_ids: int,
) -> CompartmentGrid:
    """Build per-compartment routing grid (verts, allowed mask, blocked mask)."""
    verts = ctx.compartment_hulls.get(comp)
    if verts is None or len(verts) < 3:
        x0, x1, y0, y1 = ctx.effective_region(comp)
        verts = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
    x_span = float(max(verts[:, 0].max() - verts[:, 0].min(), 1.0))
    y_span = float(max(verts[:, 1].max() - verts[:, 1].min(), 1.0))
    area = max(x_span * y_span, 1.0)
    step = float(np.clip((area / max(n_node_ids * 7.0, 64.0)) ** 0.5, 0.16, 0.45))
    pad = step * 1.5
    x_min = float(verts[:, 0].min() - pad)
    x_max = float(verts[:, 0].max() + pad)
    y_min = float(verts[:, 1].min() - pad)
    y_max = float(verts[:, 1].max() + pad)
    n_cols = int(np.ceil((x_max - x_min) / step)) + 1
    n_rows = int(np.ceil((y_max - y_min) / step)) + 1
    cols = np.arange(n_cols)
    rows = np.arange(n_rows)
    gx, gy = np.meshgrid(
        x_min + cols * step,
        y_min + rows * step,
        indexing="xy",
    )
    allowed = _polygon_mask(verts, gx, gy)
    if int(allowed.sum()) < n_node_ids:
        allowed[:] = True

    allowed_cells = np.argwhere(allowed)
    allowed_xy = np.column_stack(
        (
            x_min + allowed_cells[:, 1] * step,
            y_min + allowed_cells[:, 0] * step,
        )
    )
    blocked = ~allowed
    return CompartmentGrid(
        comp=comp,
        x_min=x_min, y_min=y_min, step=step,
        n_rows=n_rows, n_cols=n_cols,
        allowed=allowed,
        allowed_cells=allowed_cells,
        allowed_xy=allowed_xy,
        blocked=blocked,
    )


# Snap route endpoints to open cells and record routed node positions.
def _assign_node_anchors(
    grid: CompartmentGrid,
    node_ids: list[str],
    ctx: BuilderContext,
    rxn_set: set[str],
    met_set: set[str],
    routed_rxn_pos: dict[str, tuple[float, float]],
    routed_met_pos: dict[str, tuple[float, float]],
    routed_rxns: set[str],
    routed_mets: set[str],
) -> dict[str, tuple[int, int]]:
    """Snap each node to its nearest open grid cell; mutate the global routed dicts."""
    node_to_cell: dict[str, tuple[int, int]] = {}
    reserved: set[tuple[int, int]] = set()
    for node_id in node_ids:
        cell = _nearest_open_cell(
            ctx.pos[node_id], grid.allowed_cells, grid.allowed_xy, reserved,
        )
        if cell is None:
            continue
        node_to_cell[node_id] = cell
        reserved.add(cell)
        anchor_xy = grid.to_xy(cell)
        if node_id in rxn_set:
            routed_rxn_pos[node_id] = anchor_xy
            routed_rxns.add(node_id)
        if node_id in met_set:
            routed_met_pos[node_id] = anchor_xy
            routed_mets.add(node_id)
    return node_to_cell


# Route all regular edges inside a compartment grid.
def _route_compartment_edges(
    grid: CompartmentGrid,
    records: list[dict],
    node_to_cell: dict[str, tuple[int, int]],
) -> tuple[
    dict[str, list[tuple[int, int]]],
    list[str],
    dict[str, dict],
    set[tuple[int, int]],
]:
    """Run A* per record with track-reuse penalty; return route cells + bookkeeping."""
    blocked = grid.blocked
    track_usage = np.zeros((grid.n_rows, grid.n_cols), dtype=np.int16)
    route_cells: dict[str, list[tuple[int, int]]] = {}
    rec_by_lid: dict[str, dict] = {}
    ordered_lids: list[str] = []
    comp_used_cells: set[tuple[int, int]] = set()
    for rec in records:
        start = node_to_cell.get(rec["rxn_id"])
        end = node_to_cell.get(rec["met_id"])
        if start is None or end is None:
            continue
        penalty = np.where(track_usage > 0, 0.0, 0.22)
        saved = bool(blocked[start]), bool(blocked[end])
        blocked[start] = False
        blocked[end] = False
        cells = _a_star_grid(start, end, blocked, cell_penalty=penalty)
        blocked[start] = saved[0]
        blocked[end] = saved[1]
        if not cells:
            cells = [start, end]
        lid = str(rec["line_id"])
        route_cells[lid] = list(cells)
        rec_by_lid[lid] = rec
        ordered_lids.append(lid)
        for cell in cells:
            track_usage[cell] += 1
            comp_used_cells.add(cell)
    return route_cells, ordered_lids, rec_by_lid, comp_used_cells


# Generate centred lane offsets for shared rail segments.
def _lane_slot_candidates(max_lanes: int, n_routes: int) -> list[float]:
    """Route-local lane slots, ordered from visually centred outward."""
    span = max(max_lanes, n_routes, 1)
    if max_lanes <= 1:
        slots = [0.0]
        for k in range(span):
            slots.extend([-(k + 0.5), k + 0.5])
        return slots
    if max_lanes % 2 == 0:
        slots: list[float] = []
        for k in range(span):
            slots.extend([-(k + 0.5), k + 0.5])
        slots.append(0.0)
        return slots
    slots = [0.0]
    for k in range(1, span + 1):
        slots.extend([-float(k), float(k)])
    return slots


# Assign non-conflicting lane slots to routed rail paths.
def _compute_lane_slots(
    route_cells: dict[str, list[tuple[int, int]]],
    ordered_lids: list[str],
) -> dict[str, float]:
    """Pick a non-conflicting lane slot per route along its shared edges."""
    cell_edge_lines: dict[tuple, list[str]] = defaultdict(list)
    for lid in ordered_lids:
        cells = route_cells[lid]
        for i in range(len(cells) - 1):
            a, b = cells[i], cells[i + 1]
            edge = (min(a, b), max(a, b))
            if lid not in cell_edge_lines[edge]:
                cell_edge_lines[edge].append(lid)

    route_edge_sigs: dict[str, list[tuple[tuple, bool]]] = {}
    for lid in ordered_lids:
        cells = route_cells[lid]
        sigs: list[tuple[tuple, bool]] = []
        for i in range(len(cells) - 1):
            a, b = cells[i], cells[i + 1]
            edge = (min(a, b), max(a, b))
            sigs.append((edge, (a, b) == edge))
        route_edge_sigs[lid] = sigs

    used_edge_slots: dict[tuple, set[float]] = defaultdict(set)
    route_slots: dict[str, float] = {}
    for lid in ordered_lids:
        edge_sigs = route_edge_sigs[lid]
        max_lanes = max(
            (len(cell_edge_lines[edge]) for edge, _ in edge_sigs),
            default=1,
        )
        chosen_slot = 0.0
        candidates = _lane_slot_candidates(max_lanes, len(ordered_lids))
        for candidate in candidates:
            if all(
                (candidate if forward else -candidate)
                not in used_edge_slots[edge]
                for edge, forward in edge_sigs
            ):
                chosen_slot = candidate
                break
        route_slots[lid] = chosen_slot
        for edge, forward in edge_sigs:
            used_edge_slots[edge].add(
                chosen_slot if forward else -chosen_slot
            )
    return route_slots


# Serialise compartment route paths and lane slots for the frontend.
def _emit_compartment_routes(
    grid: CompartmentGrid,
    route_cells: dict[str, list[tuple[int, int]]],
    ordered_lids: list[str],
    rec_by_lid: dict[str, dict],
    route_slots: dict[str, float],
    comp: str,
    node_cells: set[tuple[int, int]],
) -> list[dict]:
    """Simplify each route's cell path into polyline + per-segment lane slots."""
    routes: list[dict] = []
    for lid in ordered_lids:
        rec = rec_by_lid[lid]
        cells = route_cells[lid]
        if len(cells) < 2:
            path_pts = [grid.to_xy(cells[0])] if cells else []
            seg_slots: list[float] = []
        else:
            leg_breaks = {0, len(cells) - 1}
            for i in range(1, len(cells) - 1):
                if cells[i] in node_cells:
                    leg_breaks.add(i)
            leg_breaks_sorted = sorted(leg_breaks)
            keep_set = set(leg_breaks)
            for li in range(len(leg_breaks_sorted) - 1):
                i0 = leg_breaks_sorted[li]
                i1 = leg_breaks_sorted[li + 1]
                for i in range(i0 + 1, i1):
                    ax, ay = cells[i - 1]
                    bx, by = cells[i]
                    cx, cy = cells[i + 1]
                    if (bx - ax) * (cy - by) != (by - ay) * (cx - bx):
                        keep_set.add(i)
            keep_idx = sorted(keep_set)
            path_pts = [grid.to_xy(cells[i]) for i in keep_idx]
            route_slot = route_slots.get(lid, 0.0)
            seg_slots = [route_slot] * (len(keep_idx) - 1)
        routes.append({
            "line_id": rec["line_id"],
            "rxn_id": rec["rxn_id"],
            "met_id": rec["met_id"],
            "comp": comp,
            "stoich": float(rec["stoich"]),
            "path": [[float(p[0]), float(p[1])] for p in path_pts],
            "seg_slots": seg_slots,
        })
    return routes


# Build dashed guide links from detached nodes to nearby tracks.
def _compute_guide_links(
    ctx: BuilderContext,
    rxn_nodes: list[str],
    met_nodes: list[str],
    comp_track_xy: dict[str, np.ndarray],
    routed_rxns: set[str],
    routed_mets: set[str],
    rxn_set: set[str],
) -> list[dict]:
    """Build dashed guide links from un-routed nodes to nearest used track cell."""
    G = ctx.G
    pos = ctx.pos
    guides: list[dict] = []
    seen: set[tuple[str, str]] = set()
    for node_id in rxn_nodes + met_nodes:
        comp = G.nodes[node_id].get("comp")
        track_xy = comp_track_xy.get(comp)
        if track_xy is None or track_xy.size == 0:
            continue
        if node_id in routed_rxns or node_id in routed_mets:
            continue
        src = pos[node_id]
        deltas = track_xy - np.asarray(src, dtype=float)
        idx = int(np.argmin(np.einsum("ij,ij->i", deltas, deltas)))
        gx_, gy_ = track_xy[idx]
        node_type = "rxn" if node_id in rxn_set else "met"
        key = (node_type, node_id)
        if key in seen:
            continue
        seen.add(key)
        guides.append({
            "node_id": node_id,
            "node_type": node_type,
            "x": [float(src[0]), float(gx_)],
            "y": [float(src[1]), float(gy_)],
        })
    return guides


# Route all regular in-compartment reaction-metabolite links.
def compute_regular_rails(
    ctx: BuilderContext,
    rxn_nodes: list[str],
    met_nodes: list[str],
    views: list[ViewSpec],
) -> dict:
    """Route regular in-compartment reaction↔metabolite links on a grid."""
    if ctx.G is None or not ctx.pos:
        return _default_rail_data(ctx, rxn_nodes, met_nodes)

    flux_priority = _compute_flux_priority(views)
    edge_records, node_degree = _collect_edge_records(
        ctx, rxn_nodes, met_nodes, flux_priority,
    )
    if not edge_records:
        return _default_rail_data(ctx, rxn_nodes, met_nodes)

    rxn_set = set(rxn_nodes)
    met_set = set(met_nodes)
    routed_rxn_pos: dict[str, tuple[float, float]] = {}
    routed_met_pos: dict[str, tuple[float, float]] = {}
    routed_rxns: set[str] = set()
    routed_mets: set[str] = set()
    routes_geom: list[dict] = []
    step_values: list[float] = []
    comp_track_xy: dict[str, np.ndarray] = {}

    for comp, records in edge_records.items():
        node_ids = sorted(
            {rec["rxn_id"] for rec in records} | {rec["met_id"] for rec in records},
            key=lambda nid: (-node_degree.get(nid, 0), nid),
        )
        if not node_ids:
            continue
        grid = _build_compartment_grid(ctx, comp, len(node_ids))
        step_values.append(grid.step)
        node_to_cell = _assign_node_anchors(
            grid, node_ids, ctx, rxn_set, met_set,
            routed_rxn_pos, routed_met_pos, routed_rxns, routed_mets,
        )

        records.sort(
            key=lambda rec: (
                -float(rec["priority"]),
                -node_degree.get(rec["rxn_id"], 0),
                -node_degree.get(rec["met_id"], 0),
                str(rec["line_id"]),
            )
        )
        comp_start_idx = len(routes_geom)
        route_cells, ordered_lids, rec_by_lid, comp_used_cells = (
            _route_compartment_edges(grid, records, node_to_cell)
        )
        route_slots = _compute_lane_slots(route_cells, ordered_lids)
        node_cells = set(node_to_cell.values())
        routes_geom.extend(_emit_compartment_routes(
            grid, route_cells, ordered_lids, rec_by_lid,
            route_slots, comp, node_cells,
        ))

        comp_track_xy[comp] = np.asarray(
            [grid.to_xy(cell) for cell in sorted(comp_used_cells)],
            dtype=float,
        )
        for ci, r in enumerate(routes_geom[comp_start_idx:]):
            r["base_color"] = _RAIL_ROUTE_PALETTE[ci % len(_RAIL_ROUTE_PALETTE)]

    for rxn_id in rxn_nodes:
        routed_rxn_pos.setdefault(rxn_id, ctx.pos[rxn_id])
    for met_id in met_nodes:
        routed_met_pos.setdefault(met_id, ctx.pos[met_id])

    guides = _compute_guide_links(
        ctx, rxn_nodes, met_nodes, comp_track_xy,
        routed_rxns, routed_mets, rxn_set,
    )

    return {
        "routes": routes_geom,
        "guides": guides,
        "rxn_pos": routed_rxn_pos,
        "met_pos": routed_met_pos,
        "routed_rxns": routed_rxns,
        "routed_mets": routed_mets,
        "grid_step": min(step_values) if step_values else 0.4,
    }


# ==========================================================================
# FIGURE HELPERS
# ==========================================================================

# Identify reactions and compartments with usable flux data.
def _compute_active_sets(
    ctx: BuilderContext,
) -> tuple[set[str], set[str]]:
    """Compute the active reaction set and the compartments they touch."""
    G = ctx.G
    cfg = ctx.cfg
    compartments = cfg["compartments"]
    exchange_key = cfg["exchange_key"]
    mdl = ctx.cobra_model

    active_rxn_set: set[str] = set()
    for view in ctx.views:
        for rxn_id, flux in view.flux_dict.items():
            if _has_flux_data(float(flux)):
                active_rxn_set.add(rxn_id)
    if not active_rxn_set:
        active_rxn_set = {
            r for r in ctx.rxn_nodes
            if ctx.rxn_kind_map.get(r) != "regular"
        }

    active_compartments: set[str] = set()
    for rxn_id in active_rxn_set:
        if rxn_id in G:
            comp = G.nodes[rxn_id].get("comp")
            if comp in compartments:
                active_compartments.add(comp)

    # Non-regular reactions render as stations between compartments —
    # include both endpoints so a compartment with only active transport
    # flux does not lose its hull.
    for rxn in mdl.reactions:
        if rxn.id not in active_rxn_set:
            continue
        kind = ctx.rxn_kind_map.get(rxn.id, "regular")
        if kind not in ("transport", "import", "export", "biomass"):
            continue
        non_exch = {_met_comp(m, cfg) for m in rxn.metabolites} - {exchange_key}
        pk = _pair_key(kind, non_exch, cfg)
        if pk is None:
            continue
        active_compartments.update(c for c in pk if c in compartments)

    return active_rxn_set, active_compartments


# Remove graph nodes that have no active flux data.
def _filter_nodes_for_drop_no_data(
    ctx: BuilderContext,
    rxn_nodes: list[str],
    met_nodes: list[str],
    met_to_rxns: dict[str, list[tuple[str, float]]],
    active_rxn_set: set[str],
    active_compartments: set[str],
) -> tuple[list[str], list[str], set[str]]:
    """Drop reactions/metabolites without flux data when ``drop_no_data`` is on."""
    G = ctx.G
    compartments = ctx.cfg["compartments"]

    rxn_nodes = [r for r in rxn_nodes if r in active_rxn_set]
    regular_remaining = set(rxn_nodes)
    active_station_mets: set[str] = set()
    for rxn_id in active_rxn_set:
        if ctx.rxn_kind_map.get(rxn_id, "regular") == "regular":
            continue
        if rxn_id not in G:
            continue
        for nb in G.neighbors(rxn_id):
            if G.nodes[nb].get("ntype") != "met":
                continue
            active_station_mets.add(nb)
            nb_comp = G.nodes[nb].get("comp")
            if nb_comp in compartments:
                active_compartments.add(nb_comp)
    met_nodes = [
        m for m in met_nodes
        if any(r in regular_remaining for r, _ in met_to_rxns.get(m, []))
        or m in active_station_mets
    ]
    return rxn_nodes, met_nodes, active_compartments


# Compute station route width limits from density scale.
def _compute_route_width_bounds(ds: float) -> tuple[float, float]:
    wmin = max(1.5, 2.0 * ds)
    wmax = max(wmin + 1.5, 7.0 * ds)
    return wmin, wmax


# Compute regular rail width limits from density scale.
def _compute_rail_width_bounds(ds: float) -> tuple[float, float]:
    wmin = max(0.9, 1.4 * ds)
    wmax = max(wmin + 0.7, 4.4 * ds)
    return wmin, wmax


# Compute global flux extents for node and route scaling.
def _compute_flux_extents(
    views: list[ViewSpec],
    routes: list[dict],
    stoich_flux: bool,
) -> tuple[float, float]:
    """Global ``abs_max_flux`` (across views) and ``abs_max_route_flux`` (rails)."""
    all_valid: list[float] = [
        float(f)
        for view in views
        for f in view.flux_dict.values()
        if _has_flux_data(float(f))
    ]
    abs_max_flux = max(
        float(np.nanmax(np.abs(all_valid))) if all_valid else 1.0, 1e-9,
    )
    all_route_valid: list[float] = []
    for view in views:
        for route in routes:
            f = view.flux_dict.get(route["rxn_id"], np.nan)
            if not _has_flux_data(float(f)):
                continue
            factor = abs(float(route["stoich"])) if stoich_flux else 1.0
            all_route_valid.append(abs(float(f)) * factor)
    abs_max_route_flux = max(
        float(np.nanmax(all_route_valid)) if all_route_valid else 1.0, 1e-9,
    )
    return abs_max_flux, abs_max_route_flux


# Return the shorter visual dimension of a compartment hull.
def _hull_short_dim(ctx: BuilderContext, comp_key: str) -> float:
    v = ctx.compartment_hulls.get(comp_key)
    if v is not None and len(v) >= 3:
        return float(min(
            v[:, 0].max() - v[:, 0].min(),
            v[:, 1].max() - v[:, 1].min(),
        ))
    return 4.0


# Build per-view reaction flux, standard deviation, and marker arrays.
def _view_per_rxn_arrays(
    view: ViewSpec, rxn_nodes: list[str], ds: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[float]]:
    fluxes = np.array(
        [view.flux_dict.get(r, np.nan) for r in rxn_nodes], dtype=float,
    )
    stds = np.array(
        [view.std_dict.get(r, np.nan) for r in rxn_nodes], dtype=float,
    )
    valid = np.array([_has_flux_data(float(f)) for f in fluxes], dtype=bool)
    log_abs = np.log1p(np.where(valid, np.abs(fluxes), 0.0))
    abs_max_local = max(
        float(np.nanmax(np.abs(fluxes))) if valid.any() else 1.0, 1e-9,
    )
    sizes = (
        log_abs / np.log1p(abs_max_local) * (26 * ds)
        + max(4, round(7 * ds))
    ).tolist()
    return fluxes, stds, valid, sizes


# Serialise reaction nodes for one frontend view.
def _build_rxn_node_list(
    ctx: BuilderContext,
    view: ViewSpec,
    rxn_nodes: list[str],
    rxn_render_pos: dict[str, tuple[float, float]],
    abs_max_flux: float,
    ds: float,
    exch_cfg: dict,
) -> list[dict]:
    G = ctx.G
    mdl = ctx.cobra_model
    compartments = ctx.cfg["compartments"]
    fluxes, stds, valid, sizes = _view_per_rxn_arrays(view, rxn_nodes, ds)

    out: list[dict] = []
    for i, r in enumerate(rxn_nodes):
        rxn_o = mdl.reactions.get_by_id(r)
        comp = G.nodes[r]["comp"]
        comp_cfg = compartments.get(comp, exch_cfg)
        comp_lbl = comp_cfg["label"]
        border_color = comp_cfg["color"]
        f, sd = float(fluxes[i]), float(stds[i])
        is_valid = bool(valid[i])
        fill_color = _flux_to_hex(f, abs_max_flux) if is_valid else "#6e7681"
        f_str = f"{f:+.4f} mmol/gDW/h" if is_valid else "--- (no data)"
        sd_str = f"{sd:.4f}" if not np.isnan(sd) else ""
        display_name = rxn_o.name or r
        subs = ", ".join(ctx.rxn_substrates.get(r, [])[:5]) or "---"
        prds = ", ".join(ctx.rxn_products.get(r, [])[:5]) or "---"
        out.append({
            "id": r,
            "name": display_name,
            "x": rxn_render_pos[r][0],
            "y": rxn_render_pos[r][1],
            "comp": comp,
            "kind": "regular",
            "flux": f if is_valid else None,
            "fill_color": fill_color,
            "border_color": border_color,
            "border_width": 1.8,
            "opacity": 0.22 if not is_valid else 0.90,
            "r": float(sizes[i]),
            "symbol": "circle",
            "hover": {
                "display_name": display_name,
                "id": r if rxn_o.name else "",
                "kind_badge": "",
                "comp_label": comp_lbl,
                "pipeline_tag": view.pipeline_tag,
                "flux_str": f_str,
                "std_str": sd_str,
                "substrates": subs,
                "products": prds,
                "extra": view.hover_extra.get(r, ""),
            },
        })
    return out


# Serialise station styling for one frontend view.
def _build_station_view_data(
    view: ViewSpec,
    stations: list[StationPair],
    abs_max_flux: float,
    route_w_bounds: tuple[float, float],
) -> list[dict]:
    wmin, wmax = route_w_bounds
    out: list[dict] = []
    for st in stations:
        line_view: list[dict] = []
        for ln in st.lines:
            flux_sum = 0.0
            for rid in ln.rxn_ids:
                f = view.flux_dict.get(rid, np.nan)
                if _has_flux_data(float(f)):
                    flux_sum += abs(float(f))
            norm = float(np.clip(flux_sum / abs_max_flux, 0.0, 1.0))
            width = wmin + (wmax - wmin) * norm
            color = (_LINE_KLASS_COLORS.get(ln.klass, "#ec4899")
                     if flux_sum > ACTIVE_EDGE_EPS
                     else _LINE_KLASS_INACTIVE)
            line_view.append({
                "klass": ln.klass,
                "flux_sum": flux_sum,
                "width": float(width),
                "color": color,
            })
        out.append({
            "pair_id": st.pair_id,
            "lines": line_view,
        })
    return out


# Serialise rail styling for one frontend view.
def _build_rail_view_data(
    view: ViewSpec,
    routes: list[dict],
    stoich_flux: bool,
    abs_max_route_flux: float,
    rail_w_bounds: tuple[float, float],
) -> list[dict]:
    _, rail_wmax = rail_w_bounds
    out: list[dict] = []
    for route in routes:
        f = view.flux_dict.get(route["rxn_id"], np.nan)
        if _has_flux_data(float(f)):
            eff_flux = float(f) * (
                abs(float(route["stoich"])) if stoich_flux else 1.0
            )
            _, mag_norm, _ = _edge_bucket(eff_flux, abs_max_route_flux)
            width = rail_wmax
            color = _modulate_route_color(
                route.get("base_color", "#58a6ff"),
                mag_norm,
                True,
            )
        else:
            width = rail_wmax
            color = "rgba(110,118,129,0.18)"
        out.append({
            "line_id": route["line_id"],
            "width": float(width),
            "color": color,
        })
    return out


# Serialise metabolite nodes for the frontend.
def _build_met_node_list(
    ctx: BuilderContext,
    met_nodes: list[str],
    met_to_rxns: dict[str, list[tuple[str, float]]],
    met_render_pos: dict[str, tuple[float, float]],
    routed_mets: set[str],
) -> list[dict]:
    G = ctx.G
    compartments = ctx.cfg["compartments"]
    out: list[dict] = []
    for n in met_nodes:
        consumers = [r for r, s in met_to_rxns.get(n, []) if s < 0]
        producers = [r for r, s in met_to_rxns.get(n, []) if s > 0]
        met_c = G.nodes[n].get("comp", "?")
        comp_cfg = compartments.get(met_c, {})
        out.append({
            "id": n,
            "name": G.nodes[n]["name"],
            "x": met_render_pos[n][0],
            "y": met_render_pos[n][1],
            "comp": met_c,
            "comp_label": comp_cfg.get("label", met_c),
            "comp_color": comp_cfg.get("color", "#8b949e"),
            "display_kind": "stop" if n in routed_mets else "detached",
            "consumers": consumers[:8],
            "producers": producers[:8],
        })
    return out


# Serialise active compartment hulls for the frontend.
def _build_compartment_list(
    ctx: BuilderContext, active_compartments: set[str],
) -> list[dict]:
    compartments = ctx.cfg["compartments"]
    out: list[dict] = []
    for comp_key, cfg_c in compartments.items():
        if ctx.drop_no_data and comp_key not in active_compartments:
            continue
        verts = ctx.compartment_hulls.get(comp_key)
        if verts is None or len(verts) < 3:
            continue
        out.append({
            "key": comp_key,
            "label": cfg_c["label"],
            "color": cfg_c["color"],
            "fill": cfg_c["fill"],
            "hull_vertices": verts.tolist(),
            "hull_tension": ctx.hull_tension,
            "label_x": float(verts[:, 0].mean()),
            "label_y": float(verts[:, 1].max()) + 0.2,
        })
    return out


# Serialise station geometry and reaction summaries for the frontend.
def _build_stations_geometry(
    ctx: BuilderContext,
    stations: list[StationPair],
    median_count: float,
    exch_cfg: dict,
) -> list[dict]:
    mdl = ctx.cobra_model
    compartments = ctx.cfg["compartments"]
    views = ctx.views
    out: list[dict] = []
    for st in stations:
        n_rxns_total = sum(len(line.rxn_ids) for line in st.lines)
        scale = float(np.log1p(n_rxns_total) / max(np.log1p(median_count), 1e-6))
        scale = float(np.clip(scale, 0.85, 1.4))
        local_size = min(
            _hull_short_dim(ctx, st.comp_a),
            _hull_short_dim(ctx, st.comp_b),
        )
        pill_h = max(0.25, 0.10 * local_size * scale)
        pill_w = 0.35 * pill_h
        n_lines = len(st.lines)
        comp_a_cfg = compartments.get(st.comp_a, exch_cfg)
        comp_b_cfg = compartments.get(st.comp_b, exch_cfg)
        line_geom: list[dict] = []
        for li, ln in enumerate(st.lines):
            lane_slot = float(li - (n_lines - 1) / 2.0)
            rxn_summaries: list[dict] = []
            for rid in ln.rxn_ids:
                rxn_o = mdl.reactions.get_by_id(rid)
                flux_per_view = []
                for view in views:
                    f = view.flux_dict.get(rid, np.nan)
                    if _has_flux_data(float(f)):
                        flux_per_view.append(f"{view.label}: {float(f):+.3f}")
                    else:
                        flux_per_view.append(f"{view.label}: ---")
                rxn_summaries.append({
                    "id": rid,
                    "name": rxn_o.name or rid,
                    "flux_per_view": flux_per_view,
                    "substrates": ctx.rxn_substrates.get(rid, [])[:5],
                    "products": ctx.rxn_products.get(rid, [])[:5],
                })
            line_geom.append({
                "klass": ln.klass,
                "rxn_ids": ln.rxn_ids,
                "rxn_summaries": rxn_summaries,
                "spine_offset_idx": ln.spine_offset_idx,
                "path": [[float(p[0]), float(p[1])] for p in st.spine],
                "seg_slots": [lane_slot] * max(len(st.spine) - 1, 0),
            })
        out.append({
            "pair_id": st.pair_id,
            "comp_a": st.comp_a, "comp_b": st.comp_b,
            "comp_a_label": comp_a_cfg["label"],
            "comp_b_label": comp_b_cfg["label"],
            "comp_a_color": comp_a_cfg["color"],
            "comp_b_color": comp_b_cfg["color"],
            "anchor_a": [st.anchor_a[0], st.anchor_a[1]],
            "anchor_b": [st.anchor_b[0], st.anchor_b[1]],
            "tangent_a": [st.tangent_a[0], st.tangent_a[1]],
            "tangent_b": [st.tangent_b[0], st.tangent_b[1]],
            "pill_height": float(pill_h),
            "pill_width": float(pill_w),
            "n_rxns": int(n_rxns_total),
            "lines": line_geom,
        })
    return out


# Compute the full x/y extent needed by frontend scales.
def _compute_data_extent(
    rxn_nodes: list[str],
    met_nodes: list[str],
    rxn_render_pos: dict[str, tuple[float, float]],
    met_render_pos: dict[str, tuple[float, float]],
    rail_data: dict,
    inner_edges: list[dict],
    stations_geometry: list[dict],
) -> tuple[list[float], list[float]]:
    all_x = (
        [rxn_render_pos[r][0] for r in rxn_nodes]
        + [met_render_pos[n][0] for n in met_nodes]
    )
    all_y = (
        [rxn_render_pos[r][1] for r in rxn_nodes]
        + [met_render_pos[n][1] for n in met_nodes]
    )
    for route in rail_data["routes"]:
        for px, py in route["path"]:
            all_x.append(float(px))
            all_y.append(float(py))
    for guide in rail_data["guides"] + inner_edges:
        all_x.extend(float(v) for v in guide["x"])
        all_y.extend(float(v) for v in guide["y"])
    for st in stations_geometry:
        all_x.extend([float(st["anchor_a"][0]), float(st["anchor_b"][0])])
        all_y.extend([float(st["anchor_a"][1]), float(st["anchor_b"][1])])
        for line in st["lines"]:
            for px, py in line["path"]:
                all_x.append(float(px))
                all_y.append(float(py))
    pad = 1.0
    x_range = [min(all_x) - pad, max(all_x) + pad] if all_x else [-12.0, 8.0]
    y_range = [min(all_y) - pad, max(all_y) + pad] if all_y else [-9.0, 7.0]
    return x_range, y_range


# Assemble the complete frontend figure payload.
def _assemble_d3_data(
    ctx: BuilderContext,
    *,
    model_name: str,
    rxn_nodes: list[str],
    met_nodes: list[str],
    rxn_render_pos: dict[str, tuple[float, float]],
    met_render_pos: dict[str, tuple[float, float]],
    rail_data: dict,
    inner_edges: list[dict],
    stations_geometry: list[dict],
    compartment_list: list[dict],
    met_node_list: list[dict],
    views_data: list[dict],
    abs_max_flux: float,
    ds: float,
    grid_step_min: float,
    x_range: list[float],
    y_range: list[float],
) -> dict:
    G = ctx.G
    mdl = ctx.cobra_model
    return {
        "meta": {
            "model_name": model_name,
            "x_range": x_range,
            "y_range": y_range,
            "height": 940,
            "abs_max_flux": float(abs_max_flux),
            "density_scale": float(ds),
            "grid_step": float(grid_step_min),
        },
        "compartments": compartment_list,
        "view_labels": [v.label for v in ctx.views],
        "views": views_data,
        "met_nodes": met_node_list,
        "rail_routes": rail_data["routes"],
        "guide_links": rail_data["guides"],
        "stations": stations_geometry,
        "inner_edges": inner_edges,
        "interactivity": {
            "view_labels": [v.label for v in ctx.views],
            "rxn_ids": rxn_nodes,
            "rxn_names": [mdl.reactions.get_by_id(r).name or r for r in rxn_nodes],
            "rxn_x": [rxn_render_pos[r][0] for r in rxn_nodes],
            "rxn_y": [rxn_render_pos[r][1] for r in rxn_nodes],
            "met_ids": met_nodes,
            "met_names": [G.nodes[n].get("name", n) for n in met_nodes],
            "met_x": [met_render_pos[n][0] for n in met_nodes],
            "met_y": [met_render_pos[n][1] for n in met_nodes],
        },
    }


# Rebuild the complete figure JSON for the widget.
def rebuild_figure(ctx: BuilderContext) -> str:
    """Build the D3 figure JSON descriptor and return it as a string."""
    cfg = ctx.cfg
    compartments = cfg["compartments"]
    exchange_key = cfg["exchange_key"]
    mdl = ctx.cobra_model

    # Drop-no-data filtering on copies. Non-regular reactions render as
    # stations; regular reactions render as in-compartment rail nodes.
    rxn_nodes = [
        r for r in ctx.rxn_nodes
        if ctx.rxn_kind_map.get(r, "regular") == "regular"
    ]
    met_nodes = list(ctx.met_nodes)
    met_to_rxns = dict(ctx.met_to_rxns)

    active_rxn_set, active_compartments = _compute_active_sets(ctx)

    if ctx.drop_no_data:
        rxn_nodes, met_nodes, active_compartments = (
            _filter_nodes_for_drop_no_data(
                ctx, rxn_nodes, met_nodes, met_to_rxns,
                active_rxn_set, active_compartments,
            )
        )

    ds = ctx.resolve_density_scale(len(rxn_nodes))
    rail_data = compute_regular_rails(ctx, rxn_nodes, met_nodes, ctx.views)
    stations, grid_step, inner_edges = compute_stations(
        ctx, active_rxn_set, active_compartments,
    )

    pair_rxn_counts = [
        sum(len(line.rxn_ids) for line in st.lines) for st in stations
    ]
    median_count = float(np.median(pair_rxn_counts)) if pair_rxn_counts else 1.0
    median_count = max(median_count, 1.0)

    route_w_bounds = _compute_route_width_bounds(ds)
    rail_w_bounds = _compute_rail_width_bounds(ds)

    abs_max_flux, abs_max_route_flux = _compute_flux_extents(
        ctx.views, rail_data["routes"], ctx.stoich_flux,
    )

    exch_cfg = compartments.get(
        exchange_key, {"color": "#8b949e", "label": exchange_key},
    )

    rxn_render_pos = rail_data["rxn_pos"]
    met_render_pos = rail_data["met_pos"]
    routed_mets = rail_data["routed_mets"]

    views_data: list[dict] = []
    for view in ctx.views:
        views_data.append({
            "label": view.label,
            "rxn_nodes": _build_rxn_node_list(
                ctx, view, rxn_nodes, rxn_render_pos, abs_max_flux, ds, exch_cfg,
            ),
            "rail_routes": _build_rail_view_data(
                view, rail_data["routes"], ctx.stoich_flux,
                abs_max_route_flux, rail_w_bounds,
            ),
            "stations": _build_station_view_data(
                view, stations, abs_max_flux, route_w_bounds,
            ),
        })

    met_node_list = _build_met_node_list(
        ctx, met_nodes, met_to_rxns, met_render_pos, routed_mets,
    )
    compartment_list = _build_compartment_list(ctx, active_compartments)
    stations_geometry = _build_stations_geometry(
        ctx, stations, median_count, exch_cfg,
    )
    x_range, y_range = _compute_data_extent(
        rxn_nodes, met_nodes, rxn_render_pos, met_render_pos,
        rail_data, inner_edges, stations_geometry,
    )

    model_name = (mdl.id or "Model") if mdl else "Model"
    d3_data = _assemble_d3_data(
        ctx,
        model_name=model_name,
        rxn_nodes=rxn_nodes,
        met_nodes=met_nodes,
        rxn_render_pos=rxn_render_pos,
        met_render_pos=met_render_pos,
        rail_data=rail_data,
        inner_edges=inner_edges,
        stations_geometry=stations_geometry,
        compartment_list=compartment_list,
        met_node_list=met_node_list,
        views_data=views_data,
        abs_max_flux=abs_max_flux,
        ds=ds,
        grid_step_min=min(grid_step, rail_data["grid_step"]),
        x_range=x_range,
        y_range=y_range,
    )

    return json.dumps(d3_data, default=_json_default)


class _BuilderMixin:
    """Compatibility adapter for the widget class.

    The builder implementation is intentionally functional: the widget takes a
    snapshot of its mutable state, then delegates to pure-ish builder functions.
    Keeping this mixin preserves the existing ``FLoV`` inheritance contract and
    older private method call sites.
    """

    # Snapshot widget state into the functional builder context.
    def _builder_context(self) -> BuilderContext:
        """Snapshot the host widget state needed by the figure builder."""
        cfg = self._cfg

        # Use the configured base region before layout has cached effective regions.
        def _fallback_effective_region(
            comp_key: str,
        ) -> tuple[float, float, float, float]:
            return cfg["compartments"][comp_key]["region"]

        return BuilderContext(
            cobra_model=self._cobra_model,
            cfg=cfg,
            G=self._G,
            pos=self._pos,
            views=self._views,
            rxn_nodes=self._rxn_nodes,
            met_nodes=self._met_nodes,
            met_to_rxns=self._met_to_rxns,
            rxn_kind_map=self._rxn_kind_map,
            rxn_substrates=self._rxn_substrates,
            rxn_products=self._rxn_products,
            compartment_hulls=getattr(self, "_compartment_hulls", {}),
            compartment_curve_samples=getattr(
                self, "_compartment_curve_samples", {},
            ),
            compartment_curve_tangents=getattr(
                self, "_compartment_curve_tangents", {},
            ),
            drop_no_data=self._drop_no_data,
            stoich_flux=self._stoich_flux,
            hull_tension=self._hull_tension,
            effective_region=getattr(
                self, "_effective_region", _fallback_effective_region,
            ),
            resolve_density_scale=self._resolve_density_scale,
        )

    # Preserve the old station-routing private method.
    def _compute_stations(
        self,
        active_rxn_set: set[str],
        active_compartments: set[str],
    ) -> tuple[list[StationPair], float, list[dict]]:
        """Delegate station routing to the functional builder."""
        return compute_stations(
            self._builder_context(), active_rxn_set, active_compartments,
        )

    # Preserve the old regular-rail private method.
    def _compute_regular_rails(
        self,
        rxn_nodes: list[str],
        met_nodes: list[str],
        views: list[ViewSpec],
    ) -> dict:
        """Delegate regular edge routing to the functional builder."""
        return compute_regular_rails(
            self._builder_context(), rxn_nodes, met_nodes, views,
        )

    # Preserve the old figure rebuild private method.
    def _rebuild_figure(self) -> None:
        """Build and store the JSON descriptor consumed by the frontend."""
        self._figure_json = rebuild_figure(self._builder_context())
