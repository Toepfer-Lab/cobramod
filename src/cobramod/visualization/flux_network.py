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
"""

from __future__ import annotations

import hashlib
import json
import re
from contextlib import contextmanager
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Union, Optional

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib  # type: ignore[no-redef]
    except ImportError:
        tomllib = None  # type: ignore[assignment]

import anywidget
import cobra
import networkx as nx
import numpy as np
import pandas as pd
from cobra import Solution
from scipy.spatial import ConvexHull
from scipy.sparse import csr_matrix
from traitlets import traitlets


# ==========================================================================
# CONSTANTS
# ==========================================================================

ACTIVE_EDGE_EPS = 1e-6
MIN_MET_RXN_DIST = 0.5
MIN_MET_EDGE_DIST = 0.5
EDGE_BINS = 6
_MET_COLOR = "#58a6ff"
EDGE_LOG_STRETCH = 99.0
EDGE_MIN_VIS_NORM = 0.22


# ==========================================================================
# DATA CLASSES
# ==========================================================================

@dataclass
class ViewSpec:
    """One flux view shown in the visualisation."""
    label: str
    flux_dict: dict[str, float]
    std_dict: dict[str, float] = field(default_factory=dict)
    hover_extra: dict[str, str] = field(default_factory=dict)
    pipeline_tag: str = ""
# ==========================================================================
# DEFAULTS  (AraCore layout — overridden by TOML config)
# ==========================================================================

_DEFAULT_COMPARTMENTS: dict[str, dict] = {
    "h": dict(region=(-11.5, -2.5, -4.5, 6.0), color="#3fb950",
              fill="rgba(63,185,80,0.15)", label="Chloroplast"),
    "c": dict(region=(-2.0, 2.0, -4.5, 6.0), color="#58a6ff",
              fill="rgba(88,166,255,0.15)", label="Cytosol"),
    "m": dict(region=(2.5, 7.5, 0.0, 6.0), color="#ffa657",
              fill="rgba(255,166,87,0.15)", label="Mitochondria"),
    "p": dict(region=(2.5, 6.5, -4.5, -0.5), color="#d2a8ff",
              fill="rgba(210,168,255,0.15)", label="Peroxisome"),
    "u": dict(region=(-11.5, 7.5, -8.5, -5.0), color="#8b949e",
              fill="rgba(139,148,158,0.10)", label="Exchange / Transport"),
}

_DEFAULT_CURRENCY_KW: set[str] = {
    "atp", "adp", "amp", "nad", "nadh", "nadp", "nadph",
    "water", "h2o", "proton", "h+", "orthophosphate", "pyrophosphate",
    "phosphate", "co2", "carbon dioxide", "oxygen", " o2",
    "coa", "coenzyme a", "fad", "fadh2",
}

_DEFAULT_IMPORT_PREFIX: tuple[str, ...] = ("Im_",)
_DEFAULT_EXPORT_PREFIX: tuple[str, ...] = ("Ex_", "Si_")
_DEFAULT_BIOMASS_PREFIX: tuple[str, ...] = ("Bio_",)
_DEFAULT_EXCHANGE_KEY = "u"


# ==========================================================================
# CONFIG LOADING
# ==========================================================================

def _load_toml(path: Path) -> dict:
    if tomllib is None:
        raise ImportError(
            "TOML config requires Python >= 3.11 or 'pip install tomli'."
        )
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def _parse_config(cfg: dict) -> dict:
    """Parse a TOML config dict into internal config format."""
    compartments: dict[str, dict] = {}
    comp_raw = cfg.get("compartment", {})
    if comp_raw:
        for key, entry in comp_raw.items():
            x0, x1, y0, y1 = [float(v) for v in entry["region"]]
            compartments[key] = dict(
                region=(x0, x1, y0, y1),
                color=entry.get("color", "#888888"),
                fill=entry.get("fill", "rgba(128,128,128,0.14)"),
                label=entry.get("label", key),
            )
    else:
        compartments = dict(_DEFAULT_COMPARTMENTS)

    lay = cfg.get("layout", {})
    exchange_key = lay.get("exchange_key", _DEFAULT_EXCHANGE_KEY)
    currency_kw = set(lay.get("currency_keywords", list(_DEFAULT_CURRENCY_KW)))
    import_prefix = tuple(lay.get("import_prefixes", list(_DEFAULT_IMPORT_PREFIX)))
    export_prefix = tuple(lay.get("export_prefixes", list(_DEFAULT_EXPORT_PREFIX)))
    biomass_prefix = tuple(lay.get("biomass_prefixes", list(_DEFAULT_BIOMASS_PREFIX)))

    # Build compartment-aware regex patterns
    _alts = "|".join(
        re.escape(k) for k in sorted(compartments.keys(), key=len, reverse=True)
    )
    suffix_re = re.compile(rf"_({_alts})(?:_[a-z])?$", re.IGNORECASE)
    bracket_re = re.compile(rf"\[({_alts})\]$", re.IGNORECASE)
    # Fallback: compartment key directly at the end of the ID without a separator
    # (e.g. "ENOc", "PDHm", "SUCOAS1m"). Sorted longest-first to avoid "m" matching
    # inside "lm" etc.
    nosep_re = re.compile(rf"({_alts})$", re.IGNORECASE)

    return dict(
        compartments=compartments,
        exchange_key=exchange_key,
        currency_kw=currency_kw,
        import_prefix=import_prefix,
        export_prefix=export_prefix,
        biomass_prefix=biomass_prefix,
        suffix_re=suffix_re,
        bracket_re=bracket_re,
        nosep_re=nosep_re,
    )


# ==========================================================================
# HELPER FUNCTIONS (stateless — take config as parameter)
# ==========================================================================

def _rxn_comp(rxn_id: str, cfg: dict) -> str:
    m = cfg["suffix_re"].search(rxn_id)
    if m:
        return m.group(1).lower()
    m = cfg["bracket_re"].search(rxn_id)
    if m:
        return m.group(1).lower()
    # Last resort: compartment key appended directly without separator (e.g. "ENOc")
    m = cfg["nosep_re"].search(rxn_id)
    return m.group(1).lower() if m else cfg["exchange_key"]

def _met_comp_str(met_id: str, cfg: dict) -> str:
    m = cfg["bracket_re"].search(met_id)
    if m:
        return m.group(1).lower()
    m = cfg["suffix_re"].search(met_id)
    return m.group(1).lower() if m else cfg["exchange_key"]


def _met_comp(met: cobra.Metabolite, cfg: dict) -> str:
    if met.compartment:
        key = met.compartment.lower()
        if key in cfg["compartments"]:
            return key
    return _met_comp_str(met.id, cfg)


def _is_currency(met: cobra.Metabolite, cfg: dict) -> bool:
    t = (met.name or "").lower() + " " + met.id.lower()
    return any(kw in t for kw in cfg["currency_kw"])


def _normalise_met_token(raw: str, cfg: dict) -> str:
    """Strip compartment suffix/bracket, collapse repeated underscores, and
    remove any residual trailing token.  Used by _is_proton on both the
    metabolite id and its name (TAIR10 stores names as 'H_plus__c0')."""
    s = cfg["bracket_re"].sub("", raw).lower().strip()
    s = cfg["suffix_re"].sub("", s).strip("_")
    s = re.sub(r"_+", "_", s)           # collapse __  →  _
    s = re.sub(r"_[a-z0-9]+$", "", s)   # strip residual compartment token
    if s.startswith("m_"):
        s = s[2:]
    return s


def _is_proton(met: cobra.Metabolite, cfg: dict) -> bool:
    name_lower = (met.name or "").lower().strip()
    # TAIR10 uses KEGG IDs (cpd00067_c0) and encodes the human name as
    # 'H_plus__c0', so we normalise *both* id and name the same way.
    id_base   = _normalise_met_token(met.id, cfg)
    name_base = _normalise_met_token(name_lower, cfg)
    _proton_tokens = {"h", "h_plus", "hplus"}
    return (
        id_base   in _proton_tokens
        or name_base in _proton_tokens
        or name_lower in {"h+", "proton", "h"}
        or "proton" in name_lower
    )


def _rxn_kind(rxn: cobra.Reaction, cfg: dict) -> str:
    if any(rxn.id.startswith(p) for p in cfg["import_prefix"]):
        return "import"
    if any(rxn.id.startswith(p) for p in cfg["export_prefix"]):
        return "export"
    if any(rxn.id.startswith(p) for p in cfg["biomass_prefix"]):
        return "biomass"
    non_exch = {_met_comp(m, cfg) for m in rxn.metabolites} - {cfg["exchange_key"]}
    if len(non_exch) > 1:
        return "transport"
    return "regular"


def _stable_unit(met_id: str) -> float:
    digest = hashlib.md5(met_id.encode("utf-8")).hexdigest()
    v = int(digest[:8], 16) / float(0xFFFFFFFF)
    return 2.0 * v - 1.0


def _has_flux_data(flux: float) -> bool:
    return (not np.isnan(flux)) and (abs(float(flux)) > ACTIVE_EDGE_EPS)


def _emit_progress(message: str) -> None:
    print(f"FLoV: {message}", flush=True)

# ==========================================================================
# LAYOUT HELPERS
# ==========================================================================

def _push_away_from_reactions(
    mx: float, my: float, rxn_nbrs: list[str], met_id: str,
    pos: dict, collision_modifier: float,
) -> tuple[float, float]:
    dist_req = MIN_MET_RXN_DIST * collision_modifier
    for rxn_id in rxn_nbrs:
        rx, ry = pos[rxn_id]
        if float(np.hypot(mx - rx, my - ry)) < dist_req:
            angle = np.pi * _stable_unit(met_id)
            return (
                rx + dist_req * float(np.cos(angle)),
                ry + dist_req * float(np.sin(angle)),
            )
    return mx, my


def _batch_point_to_segment_dist(
    px: float, py: float,
    ax: np.ndarray, ay: np.ndarray,
    bx: np.ndarray, by: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Vectorised: compute distance from point (px, py) to each of S segments.

    Parameters
    ----------
    ax, ay : (S,) arrays — segment start points
    bx, by : (S,) arrays — segment end points

    Returns
    -------
    dist, cx, cy : (S,) arrays — distance and closest point on each segment
    """
    dx, dy = bx - ax, by - ay
    seg_len_sq = dx * dx + dy * dy
    safe = seg_len_sq > 1e-12
    t = np.where(
        safe,
        np.clip(
            ((px - ax) * dx + (py - ay) * dy) / np.where(safe, seg_len_sq, 1.0),
            0.0, 1.0,
        ),
        0.0,
    )
    cx, cy = ax + t * dx, ay + t * dy
    return np.hypot(px - cx, py - cy), cx, cy


def _resolve_edge_collisions(
    met_pos: dict[str, tuple[float, float]],
    met_nodes: list[str],
    G: nx.Graph,
    pos: dict,
    collision_modifier: float,
    iterations: int = 3,
) -> dict[str, tuple[float, float]]:
    """Push metabolites away from edges that belong to other metabolites.

    The inner per-segment loop is replaced by a single call to
    `_batch_point_to_segment_dist` per metabolite per iteration.
    """
    met_rxn_adj: dict[str, set[str]] = {m: set() for m in met_nodes}
    for n1, n2 in G.edges():
        if G.nodes[n1]["ntype"] == "rxn" and G.nodes[n2]["ntype"] == "met":
            met_rxn_adj[n2].add(n1)
        elif G.nodes[n2]["ntype"] == "rxn" and G.nodes[n1]["ntype"] == "met":
            met_rxn_adj[n1].add(n2)

    # Pre-build flat segment arrays.  Segment k = edge from seg_rxn[k] → seg_met[k].
    # seg_owner[k] is the metabolite that "owns" the edge (i.e. met_pos endpoint).
    seg_owner_list: list[str] = []
    seg_rxn_list: list[str] = []
    seg_ax_list: list[float] = []
    seg_ay_list: list[float] = []

    for m2_id in met_nodes:
        for rxn_id in met_rxn_adj[m2_id]:
            if rxn_id not in pos:
                continue
            seg_owner_list.append(m2_id)
            seg_rxn_list.append(rxn_id)
            seg_ax_list.append(pos[rxn_id][0])
            seg_ay_list.append(pos[rxn_id][1])

    if not seg_owner_list:
        return met_pos

    # Use integer index arrays instead of string object arrays for fast comparisons
    met_to_int: dict[str, int] = {m: i for i, m in enumerate(met_nodes)}
    # Build a global reaction ID → integer mapping over all segments
    all_seg_rxns = list(dict.fromkeys(seg_rxn_list))  # unique, order-preserving
    rxn_to_int: dict[str, int] = {r: i for i, r in enumerate(all_seg_rxns)}

    seg_ax = np.array(seg_ax_list, dtype=float)
    seg_ay = np.array(seg_ay_list, dtype=float)
    seg_owner_int = np.array([met_to_int[m] for m in seg_owner_list], dtype=np.intp)
    seg_rxn_int = np.array([rxn_to_int[r] for r in seg_rxn_list], dtype=np.intp)

    # B endpoints (metabolite positions) — updated via fancy-index each iteration
    # met_pos_arr[i] = (x, y) of met_nodes[i]; kept in sync with met_pos dict
    met_pos_arr = np.array([met_pos[m] for m in met_nodes], dtype=float)  # (M, 2)
    seg_bx = met_pos_arr[seg_owner_int, 0].copy()
    seg_by = met_pos_arr[seg_owner_int, 1].copy()

    dist_req = MIN_MET_EDGE_DIST * collision_modifier

    for _ in range(iterations):
        # Refresh B endpoints in one fancy-index read (no Python loop)
        seg_bx[:] = met_pos_arr[seg_owner_int, 0]
        seg_by[:] = met_pos_arr[seg_owner_int, 1]

        next_arr = met_pos_arr.copy()
        for i, met_id in enumerate(met_nodes):
            mx, my = float(met_pos_arr[i, 0]), float(met_pos_arr[i, 1])
            own_rxns_int = np.array([rxn_to_int[r] for r in met_rxn_adj[met_id]
                                     if r in rxn_to_int], dtype=np.intp)

            # Mask: exclude segments owned by this metabolite or sharing its reactions
            mask = (seg_owner_int != i)
            if own_rxns_int.size:
                mask &= ~np.isin(seg_rxn_int, own_rxns_int)
            if not mask.any():
                continue

            dists, cx_arr, cy_arr = _batch_point_to_segment_dist(
                mx, my,
                seg_ax[mask], seg_ay[mask],
                seg_bx[mask], seg_by[mask],
            )

            closest = int(np.argmin(dists))
            if dists[closest] >= dist_req:
                continue

            best_cx, best_cy = float(cx_arr[closest]), float(cy_arr[closest])
            dx, dy = mx - best_cx, my - best_cy
            d = float(np.hypot(dx, dy))
            if d < 1e-9:
                angle = np.pi * _stable_unit(met_id)
                dx, dy = float(np.cos(angle)), float(np.sin(angle))
            else:
                dx /= d
                dy /= d
            next_arr[i, 0] = best_cx + dist_req * dx
            next_arr[i, 1] = best_cy + dist_req * dy

        met_pos_arr = next_arr
    # Convert back to dict
    return {met_nodes[i]: (float(met_pos_arr[i, 0]), float(met_pos_arr[i, 1]))
            for i in range(len(met_nodes))}


def _resolve_rxn_overlaps(
    pos: dict,
    rxn_nodes: list[str],
    x0: float, x1: float, y0: float, y1: float,
    min_dist: float = 0.1,
    iterations: int = 30,
    margin: float = 0.05,
) -> dict:
    """Push reaction nodes apart so they don't visually overlap, then clip
    them back inside their compartment region.

    Uses a Jacobi-style vectorised update: all pairwise displacements are
    computed from the start-of-iteration positions and applied simultaneously,
    eliminating the O(R²) Python loop.
    """
    pos = dict(pos)
    rxns = [r for r in rxn_nodes if r in pos]
    n = len(rxns)
    if n < 2:
        return pos

    coords = np.array([pos[r] for r in rxns], dtype=float)  # (N, 2)

    for _ in range(iterations):
        # diff[i, j] = coords[j] - coords[i]  →  (N, N, 2)
        diff = coords[np.newaxis, :, :] - coords[:, np.newaxis, :]
        dist = np.hypot(diff[..., 0], diff[..., 1])          # (N, N)

        needs = np.triu((dist < min_dist) & (dist > 1e-9), k=1)
        if not needs.any():
            break

        with np.errstate(divide="ignore", invalid="ignore"):
            nd = np.where(
                dist[..., np.newaxis] > 1e-9,
                diff / dist[..., np.newaxis],
                0.0,
            )  # (N, N, 2) normalised direction i→j

        # push[i,j] = half-overlap for pairs that need separating (upper tri)
        push = np.where(needs, (min_dist - dist) * 0.5, 0.0)
        push_sym = push + push.T  # make symmetric so both i→j and j→i contribute

        # displace[i] = -Σ_j  nd[i,j] * push_sym[i,j]
        # (i is pushed away from every j that is too close)
        coords -= (nd * push_sym[..., np.newaxis]).sum(axis=1)

    coords[:, 0] = np.clip(coords[:, 0], x0 + margin, x1 - margin)
    coords[:, 1] = np.clip(coords[:, 1], y0 + margin, y1 - margin)
    for i, r in enumerate(rxns):
        pos[r] = (float(coords[i, 0]), float(coords[i, 1]))
    return pos


def _fruchterman_reingold(
    G: nx.Graph,
    iterations: int = 50,
    k: Optional[float] = None,
    seed: Optional[int] = None,
) -> dict[str, tuple[float, float]]:
    """Numpy-vectorised Fruchterman-Reingold spring layout.

    Replaces ``nx.spring_layout`` for large compartments.  All pairwise
    repulsion forces are computed as a single ``(N, N)`` numpy broadcast;
    attractive forces iterate only over the edge list.  Runs ~50–200× faster
    than NetworkX's pure-Python implementation for N ≥ 200.

    Parameters
    ----------
    G : nx.Graph
        The graph to lay out.
    iterations : int
        Number of FR cooling iterations (default 50 — sufficient for layout
        quality; NetworkX default is 50, we previously used 150).
    k : float, optional
        Ideal spring length.  Defaults to ``sqrt(1 / N)``.
    seed : int, optional
        Random seed for reproducibility.
    """
    nodes = list(G.nodes())
    n = len(nodes)
    if n == 0:
        return {}
    if n == 1:
        return {nodes[0]: (0.0, 0.0)}

    rng = np.random.default_rng(seed)
    pos = rng.random((n, 2)).astype(float) * 2.0 - 1.0   # uniform in [-1, 1]²

    if k is None:
        k = float(np.sqrt(1.0 / n))
    k2 = k * k

    # Build edge index arrays once (attractive forces only touch edges).
    node_idx = {v: i for i, v in enumerate(nodes)}
    if G.number_of_edges():
        u_arr = np.array([node_idx[u] for u, _ in G.edges()], dtype=np.intp)
        v_arr = np.array([node_idx[v] for _, v in G.edges()], dtype=np.intp)
    else:
        u_arr = v_arr = np.empty(0, dtype=np.intp)

    # Scale iterations inversely with n so total O(n² × iter) work stays bounded.
    # For n=200 → 50 iters; n=1000 → 10; n=5000 → 10 (minimum floor).
    effective_iters = max(10, min(iterations, max(10, 10_000 // n)))

    t = max(0.1, 0.1 * np.sqrt(n))
    dt = t / (effective_iters + 1)

    # Pre-allocate every (N, N) array once outside the loop so we pay the
    # malloc cost once instead of once-per-iteration for 400 MB+ arrays.
    delta  = np.empty((n, n, 2), dtype=float)
    dist2  = np.empty((n, n), dtype=float)
    scalar = np.empty((n, n), dtype=float)   # k2 / (dist2 * dist) per pair
    disp   = np.empty((n, 2),  dtype=float)

    for _ in range(effective_iters):
        # ---- repulsive forces: all pairs (N, N) -------------------------
        # In-place broadcast subtract avoids re-allocating delta each iteration.
        np.subtract(pos[:, np.newaxis, :], pos[np.newaxis, :, :], out=delta)
        # dist² without creating a squared-elements temp array
        np.einsum("ijk,ijk->ij", delta, delta, out=dist2)
        np.fill_diagonal(dist2, 1.0)
        dist = np.sqrt(dist2)                # (N, N) — needed for temperature clamp

        # scalar[i,j] = k² / (dist²[i,j] * dist[i,j])  →  force = scalar * delta
        np.multiply(dist2, dist, out=scalar)
        np.divide(k2, scalar, out=scalar)
        np.fill_diagonal(scalar, 0.0)

        # disp[i] = Σ_j scalar[i,j] * delta[i,j]  — no (N,N,2) intermediate
        np.einsum("ij,ijk->ik", scalar, delta, out=disp)

        # ---- attractive forces: edges only (E,) -------------------------
        if u_arr.size:
            e_delta = pos[u_arr] - pos[v_arr]                    # (E, 2)
            e_dist = np.hypot(e_delta[:, 0], e_delta[:, 1])
            e_dist = np.maximum(e_dist, 1e-10)
            att_mag = e_dist / k
            att_disp = e_delta / e_dist[:, np.newaxis] * att_mag[:, np.newaxis]
            # np.bincount is ~20× faster than np.add.at for large edge lists
            # because it runs a vectorised C loop instead of an unbuffered loop.
            disp[:, 0] -= np.bincount(u_arr, weights=att_disp[:, 0], minlength=n)
            disp[:, 0] += np.bincount(v_arr, weights=att_disp[:, 0], minlength=n)
            disp[:, 1] -= np.bincount(u_arr, weights=att_disp[:, 1], minlength=n)
            disp[:, 1] += np.bincount(v_arr, weights=att_disp[:, 1], minlength=n)

        # ---- limit displacement by temperature --------------------------
        disp_norm = np.hypot(disp[:, 0], disp[:, 1])
        disp_norm = np.maximum(disp_norm, 1e-10)
        disp *= (np.minimum(disp_norm, t) / disp_norm)[:, np.newaxis]
        pos += disp
        t -= dt

    return {nodes[i]: (float(pos[i, 0]), float(pos[i, 1])) for i in range(n)}


def _scale_to_region(
    lp: dict, x0: float, x1: float, y0: float, y1: float, pad: float = 0.1,
) -> dict:
    xs = np.array([v[0] for v in lp.values()])
    ys = np.array([v[1] for v in lp.values()])
    xr = xs.max() - xs.min() or 1.0
    yr = ys.max() - ys.min() or 1.0
    xs_n = (xs - xs.min()) / xr
    ys_n = (ys - ys.min()) / yr
    w = (x1 - x0) * (1 - 2 * pad)
    h = (y1 - y0) * (1 - 2 * pad)
    xs_s = x0 + pad * (x1 - x0) + xs_n * w
    ys_s = y0 + pad * (y1 - y0) + ys_n * h
    return {n: (float(xs_s[i]), float(ys_s[i])) for i, n in enumerate(lp)}


def _flux_weighted_met_positions(
    edge_flux: dict[str, float],
    met_nodes: list[str],
    G: nx.Graph,
    pos: dict,
    collision_modifier: float,
    cfg: Optional[dict] = None,
    effective_region=None,
) -> dict[str, tuple[float, float]]:
    compartments = cfg["compartments"] if cfg else {}
    M = len(met_nodes)

    # All reactions that have a known position — used as the column index space.
    all_rxns = [r for r in pos if r not in set(met_nodes)]
    R = len(all_rxns)
    rxn_idx = {r: i for i, r in enumerate(all_rxns)}
    rxn_coords = np.array([[pos[r][0], pos[r][1]] for r in all_rxns], dtype=float)  # (R, 2)

    # Build a sparse (M, R) weight matrix in COO format, then convert to CSR.
    # Only same-compartment neighbours receive weight; cross-compartment
    # neighbours are excluded unless no same-compartment neighbour exists.
    # For large models (e.g. Recon3D) the dense matrix would be ~630 MB;
    # the sparse representation holds only O(M × avg_connections) entries.
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    met_comp_list: list[str] = []
    for i, met_id in enumerate(met_nodes):
        nbrs = [r for r in G.neighbors(met_id)
                if G.nodes[r]["ntype"] == "rxn" and r in rxn_idx]
        met_c = G.nodes[met_id].get("comp", "")
        met_comp_list.append(met_c)
        same_comp = [r for r in nbrs if G.nodes[r].get("comp") == met_c]
        weight_nbrs = same_comp if same_comp else nbrs
        for r in weight_nbrs:
            f = edge_flux.get(r, np.nan)
            w = (abs(float(f)) + 0.01) ** 3 if _has_flux_data(float(f)) else 0.02 ** 3
            rows.append(i)
            cols.append(rxn_idx[r])
            vals.append(w)

    W = csr_matrix(
        (np.array(vals, dtype=float), (rows, cols)),
        shape=(M, R),
    )

    # Compute all flux-weighted centroids in one sparse matrix multiply.
    row_sums = np.asarray(W.sum(axis=1)).ravel()          # (M,)
    has_weight = row_sums > 0
    centroids = np.zeros((M, 2), dtype=float)
    if has_weight.any():
        centroids[has_weight] = (
            (W[has_weight] @ rxn_coords)
            / row_sums[has_weight, np.newaxis]
        )
    # Fallback: metabolites with no weighted neighbours keep their current pos.
    for i in np.where(~has_weight)[0]:
        centroids[i] = pos.get(met_nodes[i], (0.0, 0.0))

    # Clip each metabolite to its compartment bounds.
    for i, met_c in enumerate(met_comp_list):
        if effective_region and met_c in compartments:
            x0_, x1_, y0_, y1_ = effective_region(met_c)
            centroids[i, 0] = float(np.clip(centroids[i, 0], x0_ + 0.2, x1_ - 0.2))
            centroids[i, 1] = float(np.clip(centroids[i, 1], y0_ + 0.2, y1_ - 0.2))

    met_pos: dict[str, tuple[float, float]] = {
        met_nodes[i]: (float(centroids[i, 0]), float(centroids[i, 1]))
        for i in range(M)
    }

    for met_id in met_nodes:
        rxn_nbrs = [r for r in G.neighbors(met_id)
                    if G.nodes[r]["ntype"] == "rxn" and r in pos]
        if not rxn_nbrs:
            continue
        mx, my = met_pos[met_id]
        ucx = float(np.mean([pos[r][0] for r in rxn_nbrs]))
        ucy = float(np.mean([pos[r][1] for r in rxn_nbrs]))
        dist_req = 0.22 * collision_modifier
        for rxn_id in rxn_nbrs:
            rx, ry = pos[rxn_id]
            d = float(np.hypot(mx - rx, my - ry))
            if d < dist_req:
                pull_dx, pull_dy = ucx - rx, ucy - ry
                pull_d = float(np.hypot(pull_dx, pull_dy))
                if pull_d > 1e-6:
                    mx = rx + dist_req * (pull_dx / pull_d)
                    my = ry + dist_req * (pull_dy / pull_d)
                else:
                    if d > 1e-9:
                        mx = rx + dist_req * ((mx - rx) / d)
                        my = ry + dist_req * ((my - ry) / d)
                    else:
                        angle = np.pi * _stable_unit(met_id)
                        mx = rx + dist_req * float(np.cos(angle))
                        my = ry + dist_req * float(np.sin(angle))
        # Re-clip after reaction-collision nudge to keep metabolite in its compartment
        met_c2 = G.nodes[met_id].get("comp", "")
        if effective_region and met_c2 in compartments:
            _x0, _x1, _y0, _y1 = effective_region(met_c2)
            mx = float(np.clip(mx, _x0 + 0.2, _x1 - 0.2))
            my = float(np.clip(my, _y0 + 0.2, _y1 - 0.2))
        met_pos[met_id] = (mx, my)

    return met_pos


# ==========================================================================
# EDGE HELPERS
# ==========================================================================

def _edge_bucket(flux: float, abs_max: float) -> tuple[int, float, int]:
    abs_flux = abs(flux)
    if abs_max > 1e-9:
        rel = abs_flux / abs_max
        mag_norm = np.log1p(EDGE_LOG_STRETCH * rel) / np.log1p(EDGE_LOG_STRETCH)
        if abs_flux > 0.0:
            mag_norm = EDGE_MIN_VIS_NORM + (1.0 - EDGE_MIN_VIS_NORM) * mag_norm
    else:
        mag_norm = 0.0
    b = int(np.clip(np.floor(mag_norm * EDGE_BINS), 0, EDGE_BINS - 1))
    sign = 1 if flux >= 0 else -1
    return b, float(mag_norm), sign


def _edge_rgba(mag_norm: float, sign: int, has_flux: bool, density_scale: float = 1.0, color_modifier: float = 1.0) -> str:
    if not has_flux:
        return "rgba(110,118,129,0.20)"
    # Higher alpha on dark background — vivid colours don't need to be subdued
    alpha = (0.55 + 0.45 * mag_norm) * (0.60 + 0.40 * density_scale)
    eff_norm = float(np.clip(mag_norm * color_modifier, 0.0, 1.0))
    if sign > 0:
        # Positive: warm orange → vivid red
        r = int(255)
        g = int(180 - 180 * eff_norm)
        b_ = int(70 - 70 * eff_norm)
    else:
        # Negative: vivid sky-blue → deep blue
        r = int(100 - 50 * eff_norm)
        g = int(190 - 90 * eff_norm)
        b_ = int(255)
    return f"rgba({r},{g},{b_},{alpha:.3f})"


def _flux_to_hex(flux: float, abs_max: float) -> str:
    """Map a flux value to a dark-mode hex colour using a diverging red/blue scale."""
    if abs_max < 1e-9 or not _has_flux_data(flux):
        return "#6e7681"
    t = float(np.clip(abs(flux) / abs_max, 0.0, 1.0))
    if flux >= 0:
        # Low positive: #ef9a9a  High positive: #c62828
        r = int(239 + (198 - 239) * t)
        g = int(154 + (40 - 154) * t)
        b_ = int(154 + (40 - 154) * t)
    else:
        # Low negative: #90caf9  High negative: #1565c0
        r = int(144 + (21 - 144) * t)
        g = int(202 + (101 - 202) * t)
        b_ = int(249 + (192 - 249) * t)
    return f"#{r:02x}{g:02x}{b_:02x}"


def _edge_coords(
    pairs: list[tuple[str, str]],
    met_positions: dict[str, tuple[float, float]],
    pos: dict,
) -> tuple[list[float | None], list[float | None], list[str]]:
    ex, ey = [], []
    met_ids: list[str] = []
    for rxn_id, met_id in pairs:
        x0_, y0_ = pos[rxn_id]
        x1_, y1_ = met_positions.get(met_id, pos.get(met_id, (0.0, 0.0)))
        ex.extend([x0_, x1_, None])
        ey.extend([y0_, y1_, None])
        met_ids.append(met_id)
    return ex, ey, met_ids


def _build_edge_coords_for_view(
    edge_flux: dict[str, float],
    met_positions: dict[str, tuple[float, float]],
    oriented_edges: list[tuple[str, str, str, float]],
    pos: dict,
    scale_by_stoich: bool = False,
) -> tuple[list[list[float | None]], list[list[float | None]], list[list[str]]]:
    if scale_by_stoich:
        scaled = [abs(stoich) * abs(float(edge_flux[rxn_id]))
                  for _, rxn_id, _, stoich in oriented_edges
                  if rxn_id in edge_flux and _has_flux_data(float(edge_flux[rxn_id]))]
        abs_max = max(float(np.nanmax(scaled)) if scaled else 1.0, 1e-9)
    else:
        valid = [abs(v) for v in edge_flux.values() if _has_flux_data(float(v))]
        abs_max = max(float(np.nanmax(valid)) if valid else 1.0, 1e-9)

    buckets: dict[tuple[str, int, int], list[tuple[str, str]]] = {}
    for edge_class, rxn_id, met_id, stoich in oriented_edges:
        f = edge_flux.get(rxn_id, np.nan)
        if not _has_flux_data(float(f)):
            continue
        f_eff = float(f) * abs(stoich) if scale_by_stoich else float(f)
        b, _, sign = _edge_bucket(f_eff, abs_max)
        buckets.setdefault((edge_class, b, sign), []).append((rxn_id, met_id))

    x_lists: list[list[float | None]] = []
    y_lists: list[list[float | None]] = []
    met_id_lists: list[list[str]] = []
    for edge_class in ("reg", "ss"):
        for sign in (1, -1):
            for b in range(EDGE_BINS):
                pairs = buckets.get((edge_class, b, sign), [])
                ex, ey, met_ids = _edge_coords(pairs, met_positions, pos)
                x_lists.append(ex if ex else [None])
                y_lists.append(ey if ey else [None])
                met_id_lists.append(met_ids)
    return x_lists, y_lists, met_id_lists


# ==========================================================================
# COLOUR HELPERS
# ==========================================================================

def _tint_hex(hex_color: str, alpha: float = 0.15) -> str:
    h = hex_color.lstrip("#")
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"


# ==========================================================================
# SHAPE HELPERS
# ==========================================================================

def _smooth_hull_path(verts: np.ndarray, tension: float = 0.4) -> str:
    n = len(verts)
    parts = [f"M {verts[0][0]:.4f},{verts[0][1]:.4f}"]
    for i in range(n):
        p0 = verts[(i - 1) % n]
        p1 = verts[i % n]
        p2 = verts[(i + 1) % n]
        p3 = verts[(i + 2) % n]
        cp1 = p1 + tension * (p2 - p0) / 3.0
        cp2 = p2 - tension * (p3 - p1) / 3.0
        parts.append(
            f"C {cp1[0]:.4f},{cp1[1]:.4f} "
            f"{cp2[0]:.4f},{cp2[1]:.4f} "
            f"{p2[0]:.4f},{p2[1]:.4f}"
        )
    parts.append("Z")
    return " ".join(parts)


# ==========================================================================
# TRACE BUILDING
# ==========================================================================

def _symbol_list(
    fluxes: np.ndarray, rxn_nodes: list[str], rxn_kind_map: dict[str, str],
) -> list[str]:
    syms = []
    for i, r in enumerate(rxn_nodes):
        k = rxn_kind_map.get(r, "regular")
        if k == "import":
            syms.append("triangle-up")
        elif k in ("export", "biomass"):
            f = fluxes[i]
            syms.append("cross" if not _has_flux_data(float(f)) else "triangle-down")
        elif k == "transport":
            syms.append("hexagon2")
        else:
            syms.append("circle")
    return syms


def _make_rxn_trace(
    spec: ViewSpec,
    rxn_nodes: list[str],
    pos: dict,
    G: nx.Graph,
    cobra_model: cobra.Model,
    compartments: dict[str, dict],
    exchange_key: str,
    rxn_kind_map: dict[str, str],
    rxn_substrates: dict[str, list[str]],
    rxn_products: dict[str, list[str]],
    density_scale: float = 1.0,
) -> go.Scattergl:
    fluxes = np.array([spec.flux_dict.get(r, np.nan) for r in rxn_nodes], dtype=float)
    stds = np.array([spec.std_dict.get(r, np.nan) for r in rxn_nodes], dtype=float)
    valid = np.array([_has_flux_data(float(f)) for f in fluxes], dtype=bool)

    abs_max = max(float(np.nanmax(np.abs(fluxes))) if valid.any() else 1.0, 1e-9)
    log_abs = np.log1p(np.where(valid, np.abs(fluxes), 0.0))
    sizes = (log_abs / np.log1p(abs_max) * (26 * density_scale) + max(4, round(7 * density_scale))).tolist()

    _exch_cfg = compartments.get(exchange_key, {"color": "#888888", "label": exchange_key})
    border_colors = [
        compartments.get(G.nodes[r]["comp"], _exch_cfg)["color"]
        for r in rxn_nodes
    ]
    border_widths = [
        2.5 if rxn_kind_map.get(r, "regular") != "regular" else 1.8
        for r in rxn_nodes
    ]
    opacities = [0.22 if not v else 0.90 for v in valid]

    hover = []
    bg_colors = [_tint_hex(c) for c in border_colors]
    for i, r in enumerate(rxn_nodes):
        f, sd = fluxes[i], stds[i]
        rxn_o = cobra_model.reactions.get_by_id(r)
        comp = G.nodes[r]["comp"]
        comp_lbl = compartments.get(comp, _exch_cfg)["label"]
        kind = rxn_kind_map.get(r, "regular")

        f_str = f"{f:+.4f} mmol/gDW/h" if _has_flux_data(float(f)) else "--- (no data)"
        sd_str = f"<br><i>std = {sd:.4f}</i>" if not np.isnan(sd) else ""
        kind_badge = f" <b>[{kind.upper()}]</b>" if kind != "regular" else ""
        pipe_str = f" | Pipeline: {spec.pipeline_tag}" if spec.pipeline_tag else ""
        extra = spec.hover_extra.get(r, "")

        subs = ", ".join(rxn_substrates.get(r, [])[:5]) or "---"
        prds = ", ".join(rxn_products.get(r, [])[:5]) or "---"

        display_name = rxn_o.name or r
        id_str = f" <i>({r})</i>" if rxn_o.name else ""
        hover.append(
            f"<b>{display_name}</b>{id_str}{kind_badge}<br>"
            f"<span style='color:#888'>─────────────────────────</span><br>"
            f"<i>{comp_lbl}</i>{pipe_str}<br>"
            f"Flux: {f_str}{sd_str}<br>"
            f"Substrates: {subs}<br>"
            f"Products: {prds}"
            f"{extra}"
            f"<extra></extra>"
        )

    return go.Scattergl(
        x=[pos[r][0] for r in rxn_nodes],
        y=[pos[r][1] for r in rxn_nodes],
        mode="markers",
        marker=dict(
            size=sizes,
            symbol=_symbol_list(fluxes, rxn_nodes, rxn_kind_map),
            color=np.where(valid, fluxes, 0.0).tolist(),
            colorscale="RdBu_r",
            cmin=-abs_max, cmid=0, cmax=abs_max,
            colorbar=dict(
                title=dict(text="Flux<br>(mmol/gDW/h)",
                           side="right", font=dict(size=11)),
                x=-0.12, xanchor="left",
                thickness=14, len=0.50, y=0.5,
                tickfont=dict(size=9), outlinewidth=0,
            ),
            line=dict(color=border_colors, width=border_widths),
            opacity=opacities,
        ),
        customdata=rxn_nodes,
        hovertemplate=hover,
        hoverlabel=dict(
            bgcolor=bg_colors,
            bordercolor=border_colors,
            font=dict(color="#2c3e50", size=12),
        ),
        name=spec.label, showlegend=False,
    )


# ==========================================================================
# CSV LOADER
# ==========================================================================

def _load_ignore_csv(path: Union[str, Path]) -> tuple[set[str], set[str]]:
    """Parse a CSV listing reaction and metabolite IDs to exclude from the graph.

    Expected format — two columns, order does not matter, extra columns ignored::

        reactions,metabolites
        RXN_001,cpd00067_c0
        RXN_002,cpd00009_m0
        ,cpd00012_d0

    Empty cells are skipped.  Column headers are matched case-insensitively;
    accepted names are ``reactions``/``rxn``/``rxn_id`` and
    ``metabolites``/``met``/``met_id``.

    Returns
    -------
    rxn_ids : set[str]
        Reaction IDs to ignore.
    met_ids : set[str]
        Metabolite IDs to ignore.
    """
    df = pd.read_csv(Path(path), dtype=str)
    df.columns = df.columns.str.strip().str.lower()

    _rxn_kw = ["reactions", "rxn", "rxn_id", "reaction"]
    _met_kw = ["metabolites", "met", "met_id", "metabolite"]

    rxn_col = next((c for c in _rxn_kw if c in df.columns), None)
    met_col = next((c for c in _met_kw if c in df.columns), None)

    rxn_ids: set[str] = set()
    met_ids: set[str] = set()

    if rxn_col:
        rxn_ids = {v.strip() for v in df[rxn_col].dropna() if v.strip()}
    if met_col:
        met_ids = {v.strip() for v in df[met_col].dropna() if v.strip()}

    return rxn_ids, met_ids


def _load_single_csv(path: Path) -> pd.Series:
    df = pd.read_csv(path)
    df.columns = df.columns.str.strip()
    _id_kw = ["reaction_id", "rxn_id", "rxn", "id", "reaction"]
    _val_kw = ["flux", "value", "fba_flux", "v"]
    id_col = next((c for c in _id_kw if c in df.columns), df.columns[0])
    val_col = next((c for c in _val_kw if c in df.columns), df.columns[1])
    return pd.Series(
        pd.to_numeric(df[val_col], errors="coerce").values,
        index=df[id_col].astype(str).str.strip(),
    )


# ==========================================================================
# MAIN WIDGET CLASS
# ==========================================================================

class FLoV(anywidget.AnyWidget):
    """
    Interactive 2D Plotly-based metabolic flux network widget for Jupyter.

    Renders a bipartite graph (reactions <-> metabolites) with
    compartment-aware layout, flux-coloured edges and nodes,
    view switching, search, and metabolite dragging.

    Usage::

        widget = FLoV()
        widget.compartments = "config/aracore.toml"  # or dict
        widget.model = cobra.io.read_sbml_model("model.xml")
        widget.add_view("FBA", model.optimize())
        widget  # display in Jupyter
    """

    # anywidget trait synced to frontend
    _figure_json = traitlets.Unicode("{}").tag(sync=True)

    def __init__(
        self,
        drop_no_data: bool = False,
        drop_protons: bool = False,
        stoich_flux: bool = False,
        scale_compartments: bool = False,
        hull_tension: float = 0.4,
        collision_modifier: float = 0.001,
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
        self._stoich_flux = stoich_flux
        self._scale_compartments = scale_compartments
        self._hull_tension = hull_tension
        self._collision_modifier = collision_modifier
        self._density_scale = density_scale
        self._color_modifier = color_modifier
        self._ignored_rxns: set[str] = set()
        self._ignored_mets: set[str] = set()
        self._batch_depth: int = 0
        self._batch_needs: set[str] = set()
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
            raw = _load_toml(Path(value))
            self._cfg = _parse_config(raw)
        elif isinstance(value, dict):
            self._cfg = _parse_config(value)
        else:
            raise TypeError(f"Expected str, Path, or dict, got {type(value)}")
        self._trigger("graph")

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

    @property
    def drop_no_data(self) -> bool:
        """Remove reactions with no flux data in any view."""
        return self._drop_no_data

    @drop_no_data.setter
    def drop_no_data(self, value: bool) -> None:
        self._drop_no_data = value
        self._trigger("figure")

    @property
    def stoich_flux(self) -> bool:
        """Scale edge width by stoichiometric coefficient × flux."""
        return self._stoich_flux

    @stoich_flux.setter
    def stoich_flux(self, value: bool) -> None:
        self._stoich_flux = value
        self._trigger("figure")

    @property
    def hull_tension(self) -> float:
        """Smoothness of compartment hull outlines (0=sharp, 1=very round)."""
        return self._hull_tension

    @hull_tension.setter
    def hull_tension(self, value: float) -> None:
        self._hull_tension = value
        self._trigger("figure")

    @property
    def scale_compartments(self) -> bool:
        """Scale compartment regions proportionally to sqrt(n_reactions)."""
        return self._scale_compartments

    @scale_compartments.setter
    def scale_compartments(self, value: bool) -> None:
        self._scale_compartments = value
        self._trigger("layout")

    @property
    def collision_modifier(self) -> float:
        """Scale factor for node collision avoidance."""
        return self._collision_modifier

    @collision_modifier.setter
    def collision_modifier(self, value: float) -> None:
        self._collision_modifier = value
        self._trigger("layout")

    @property
    def density_scale(self) -> Optional[float]:
        """Visual density scale in [0, 1] applied to marker sizes, edge widths, and
        colour alpha. ``None`` (default) auto-computes the scale from network size,
        shrinking visuals for models with more than ~80 rendered reactions."""
        return self._density_scale

    @density_scale.setter
    def density_scale(self, value: Optional[float]) -> None:
        self._density_scale = value
        self._trigger("figure")

    @property
    def color_modifier(self) -> float:
        """Scale factor for edge colour saturation. 1.0 = default; >1 more vivid, <1 more faded."""
        return self._color_modifier

    @color_modifier.setter
    def color_modifier(self, value: float) -> None:
        self._color_modifier = value
        self._trigger("figure")

    @contextmanager
    def configure(self):
        """Suppress all rebuilds while setting configuration properties.

        Without this, each assignment on an already-built widget triggers its
        own expensive rebuild immediately.  Wrap configuration assignments in
        this block, then trigger one rebuild explicitly afterwards::

            with widget.configure():
                widget.drop_protons = True
                widget.density_scale = 0.6
                widget.collision_modifier = 0.5
            widget.model = model  # single rebuild here
        """
        self._batch_depth += 1
        try:
            yield
        finally:
            self._batch_depth -= 1
            # Discard queued flags — the caller is responsible for triggering
            # the rebuild explicitly (e.g. widget.model = model) after the block.
            if self._batch_depth == 0:
                self._batch_needs.clear()

    def _trigger(self, level: str) -> None:
        """Schedule or immediately execute a rebuild at the given level.

        ``level`` is one of ``'figure'``, ``'layout'``, or ``'graph'``
        (each implies all cheaper levels).  When inside a ``configure()``
        block the rebuild is deferred; otherwise it runs immediately.
        """
        if not self._graph_built:
            return
        if self._batch_depth:
            self._batch_needs.add(level)
            return
        if level == "graph":
            self._build_graph()
        elif level == "layout":
            self._compute_layout()
        self._rebuild_figure()

    def _resolve_density_scale(self, n_rxns: int) -> float:
        """Return the effective density scale: user value if set, else auto from network size."""
        if self._density_scale is not None:
            return float(np.clip(self._density_scale, 0.0, 1.0))
        return max(0.40, min(1.0, (80.0 / max(80, n_rxns)) ** 0.5))

    @property
    def drop_protons(self) -> bool:
        """Exclude proton (H+) metabolites from the graph."""
        return self._drop_protons

    @drop_protons.setter
    def drop_protons(self, value: bool) -> None:
        self._drop_protons = value
        self._trigger("graph")

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
        self._trigger("graph")

    # -- Public API --

    def add_view(
        self,
        label: str,
        flux_data: Union[Solution, dict, pd.Series, str, Path],
        std_data: Optional[dict] = None,
        pipeline_tag: str = "",
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
        )
        for i, existing in enumerate(self._views):
            if existing.label == label:
                self._views[i] = spec
                break
        else:
            self._views.append(spec)
        if self._graph_built:
            _emit_progress(f"Rendering view '{label}'...")
            self._rebuild_figure()

    def clear_views(self) -> None:
        """Remove all flux views."""
        self._views.clear()
        self._figure_json = "{}"

    def run_fba(self, label: str = "FBA") -> None:
        """Run FBA on the model and add as a view."""
        if self._cobra_model is None:
            raise ValueError("No model set.")
        _emit_progress("Running FBA...")
        sol = self._cobra_model.optimize()
        if sol.status != "optimal":
            raise RuntimeError(f"FBA status: {sol.status}")
        self.add_view(label, sol)

    def run_pfba(self, label: str = "pFBA") -> None:
        """Run parsimonious FBA on the model and add as a view."""
        if self._cobra_model is None:
            raise ValueError("No model set.")
        _emit_progress("Running pFBA...")
        sol = cobra.flux_analysis.pfba(self._cobra_model)
        self.add_view(label, sol)

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
            for met in rxn.metabolites:
                is_p = _is_proton(met, cfg)
                if not _is_currency(met, cfg) or (is_p and not self._drop_protons):
                    met_rxn_count[met.id] = met_rxn_count.get(met.id, 0) + 1

        self._met_to_rxns = {}
        for rxn in mdl.reactions:
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
            if self._rxn_kind_map[rxn.id] == "regular":
                comp = _rxn_comp(rxn.id, cfg)
                # When neither suffix_re, bracket_re, nor nosep_re could extract a
                # compartment (all return exchange_key), derive it from the majority
                # compartment of the reaction's own metabolites.  This is the most
                # reliable signal for genome-scale models like Recon3D whose reaction
                # IDs carry no compartment suffix at all.
                if comp == cfg["exchange_key"]:
                    comp_counts: dict[str, int] = {}
                    for met_r in rxn.metabolites:
                        mc_r = _met_comp(met_r, cfg)
                        comp_counts[mc_r] = comp_counts.get(mc_r, 0) + 1
                    if comp_counts:
                        # Prefer non-exchange compartments; exchange only as last resort
                        non_exch = {k: v for k, v in comp_counts.items()
                                    if k != cfg["exchange_key"]}
                        pool = non_exch if non_exch else comp_counts
                        comp = max(pool, key=pool.__getitem__)
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
            G.add_node(met_id, ntype="met", comp=mc, name=met_obj.name or met_id)
            for rxn_id, stoich in rxn_list:
                if rxn_id not in self._ignored_rxns:
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
            pos.update(_resolve_rxn_overlaps(
                pos, rxns_here, x0, x1, y0, y1,
                min_dist=0.4 * self._collision_modifier,
            ))

        # Metabolite positions — same-compartment-preferred uniform centroid,
        # computed with numpy means over batches of neighbour coordinates,
        # followed by compartment clipping.
        for i, met_id in enumerate(self._met_nodes):
            nbrs = [n for n in G.neighbors(met_id)
                    if G.nodes[n]["ntype"] == "rxn" and n in pos]
            met_c = G.nodes[met_id].get("comp", cfg["exchange_key"])
            if nbrs:
                same = [n for n in nbrs if G.nodes[n].get("comp") == met_c]
                active = same if same else nbrs
                pts = np.array([[pos[n][0], pos[n][1]] for n in active], dtype=float)
                cx, cy = float(pts[:, 0].mean()), float(pts[:, 1].mean())
                cx, cy = _push_away_from_reactions(
                    cx, cy, nbrs, met_id, pos, self._collision_modifier
                )
                pos[met_id] = (cx, cy)
            else:
                pos[met_id] = (0.0, 0.0)

            if met_c in compartments:
                x0, x1, y0, y1 = _effective_region(met_c)
                mx, my = pos[met_id]
                pos[met_id] = (
                    float(np.clip(mx, x0 + 0.2, x1 - 0.2)),
                    float(np.clip(my, y0 + 0.2, y1 - 0.2)),
                )

        # Resolve edge collisions
        _emit_progress("Resolving edge collisions...")
        _base_resolved = _resolve_edge_collisions(
            {m: pos[m] for m in self._met_nodes},
            self._met_nodes, G, pos, self._collision_modifier,
        )
        for _m in self._met_nodes:
            pos[_m] = _base_resolved[_m]
            # Re-clip: edge-collision resolver can push metabolites across compartment
            # boundaries; enforce hard bounds a second time.
            _mc = G.nodes[_m].get("comp", cfg["exchange_key"])
            if _mc in compartments:
                _x0, _x1, _y0, _y1 = _effective_region(_mc)
                _mx, _my = pos[_m]
                pos[_m] = (
                    float(np.clip(_mx, _x0 + 0.2, _x1 - 0.2)),
                    float(np.clip(_my, _y0 + 0.2, _y1 - 0.2)),
                )

        self._pos = pos

        # Store effective region function for later use
        self._effective_region = _effective_region

    # -- Internal: figure building --

    def _rebuild_figure(self) -> None:
        """Rebuild the D3 figure JSON and sync to frontend."""
        if not self._graph_built or not self._views:
            return

        _emit_progress("Rendering figure...")

        G = self._G
        pos = self._pos
        cfg = self._cfg
        compartments = cfg["compartments"]
        exchange_key = cfg["exchange_key"]
        views = self._views
        mdl = self._cobra_model

        # Drop-no-data filtering (on copies)
        rxn_nodes = list(self._rxn_nodes)
        met_nodes = list(self._met_nodes)
        met_to_rxns = dict(self._met_to_rxns)

        if self._drop_no_data:
            rxns_with_data: set[str] = set()
            for view in views:
                for rxn_id, flux in view.flux_dict.items():
                    if _has_flux_data(float(flux)):
                        rxns_with_data.add(rxn_id)
            rxn_nodes = [r for r in rxn_nodes if r in rxns_with_data]
            remaining = set(rxn_nodes)
            met_nodes = [m for m in met_nodes if any(
                r in remaining for r, _ in met_to_rxns.get(m, [])
            )]

        # Flux-weighted metabolite positions per view
        flux_met_pos_by_view = [
            _flux_weighted_met_positions(
                v.flux_dict, met_nodes, G, pos, self._collision_modifier,
                cfg=cfg, effective_region=self._effective_region,
            )
            for v in views
        ]

        # Build oriented edges
        rxn_set = set(rxn_nodes)
        met_set = set(met_nodes)
        oriented_edges: list[tuple[str, str, str, float]] = []
        for n1, n2, edata in G.edges(data=True):
            if G.nodes[n1]["ntype"] == "rxn" and G.nodes[n2]["ntype"] == "met":
                rxn_id, met_id = n1, n2
            elif G.nodes[n2]["ntype"] == "rxn" and G.nodes[n1]["ntype"] == "met":
                rxn_id, met_id = n2, n1
            else:
                continue
            if rxn_id not in pos or met_id not in pos:
                continue
            if rxn_id not in rxn_set or met_id not in met_set:
                continue
            edge_class = ("ss" if self._rxn_kind_map.get(rxn_id, "regular") != "regular"
                          else "reg")
            stoich = float(edata.get("stoich", 1.0))
            oriented_edges.append((edge_class, rxn_id, met_id, stoich))

        # Edge coordinates per view
        edge_flux_x_by_view: list[list[list[float | None]]] = []
        edge_flux_y_by_view: list[list[list[float | None]]] = []
        for vi, view in enumerate(views):
            exf, eyf, _ = _build_edge_coords_for_view(
                view.flux_dict, flux_met_pos_by_view[vi],
                oriented_edges, pos, self._stoich_flux,
            )
            edge_flux_x_by_view.append(exf)
            edge_flux_y_by_view.append(eyf)

        # Edge visual styles — same bucket ordering as _build_edge_coords_for_view
        ds = self._resolve_density_scale(len(rxn_nodes))
        _edge_wmin = max(0.5, 1.2 * ds)
        _edge_wmax = max(_edge_wmin + 0.4, 4.6 * ds)
        edge_styles: list[dict] = []
        for edge_class in ("reg", "ss"):
            for sign in (1, -1):
                for b in range(EDGE_BINS):
                    mag_norm = (b + 0.5) / EDGE_BINS
                    width = _edge_wmin + (_edge_wmax - _edge_wmin) * mag_norm
                    color = _edge_rgba(mag_norm=mag_norm, sign=sign, has_flux=True,
                                       density_scale=ds, color_modifier=self._color_modifier)
                    edge_styles.append({"class": edge_class, "color": color, "width": width})

        # Global abs_max across all views (consistent colorscale)
        all_valid: list[float] = [
            float(f)
            for view in views
            for f in view.flux_dict.values()
            if _has_flux_data(float(f))
        ]
        abs_max_flux = max(float(np.nanmax(np.abs(all_valid))) if all_valid else 1.0, 1e-9)

        n_views = len(views)
        _exch_cfg = compartments.get(exchange_key, {"color": "#8b949e", "label": exchange_key})

        # Per-view data
        views_data: list[dict] = []
        for vi, view in enumerate(views):
            fluxes = np.array([view.flux_dict.get(r, np.nan) for r in rxn_nodes], dtype=float)
            stds = np.array([view.std_dict.get(r, np.nan) for r in rxn_nodes], dtype=float)
            valid = np.array([_has_flux_data(float(f)) for f in fluxes], dtype=bool)
            log_abs = np.log1p(np.where(valid, np.abs(fluxes), 0.0))
            abs_max_local = max(float(np.nanmax(np.abs(fluxes))) if valid.any() else 1.0, 1e-9)
            sizes = (log_abs / np.log1p(abs_max_local) * (26 * ds) + max(4, round(7 * ds))).tolist()

            rxn_node_list: list[dict] = []
            for i, r in enumerate(rxn_nodes):
                rxn_o = mdl.reactions.get_by_id(r)
                comp = G.nodes[r]["comp"]
                comp_cfg = compartments.get(comp, _exch_cfg)
                comp_lbl = comp_cfg["label"]
                border_color = comp_cfg["color"]
                kind = self._rxn_kind_map.get(r, "regular")
                f, sd = float(fluxes[i]), float(stds[i])
                is_valid = bool(valid[i])
                fill_color = _flux_to_hex(f, abs_max_flux) if is_valid else "#6e7681"
                f_str = f"{f:+.4f} mmol/gDW/h" if is_valid else "--- (no data)"
                sd_str = f"{sd:.4f}" if not np.isnan(sd) else ""
                kind_badge = kind.upper() if kind != "regular" else ""
                display_name = rxn_o.name or r
                subs = ", ".join(self._rxn_substrates.get(r, [])[:5]) or "---"
                prds = ", ".join(self._rxn_products.get(r, [])[:5]) or "---"
                if kind in ("export", "biomass"):
                    symbol = "cross" if not is_valid else "triangle-down"
                elif kind == "import":
                    symbol = "triangle-up"
                elif kind == "transport":
                    symbol = "hexagon"
                else:
                    symbol = "circle"
                rxn_node_list.append({
                    "id": r,
                    "name": display_name,
                    "x": pos[r][0],
                    "y": pos[r][1],
                    "comp": comp,
                    "kind": kind,
                    "flux": f if is_valid else None,
                    "fill_color": fill_color,
                    "border_color": border_color,
                    "border_width": 2.5 if kind != "regular" else 1.8,
                    "opacity": 0.22 if not is_valid else 0.90,
                    "r": float(sizes[i]),
                    "symbol": symbol,
                    "hover": {
                        "display_name": display_name,
                        "id": r if rxn_o.name else "",
                        "kind_badge": kind_badge,
                        "comp_label": comp_lbl,
                        "pipeline_tag": view.pipeline_tag,
                        "flux_str": f_str,
                        "std_str": sd_str,
                        "substrates": subs,
                        "products": prds,
                        "extra": view.hover_extra.get(r, ""),
                    },
                })

            edge_groups: list[dict] = []
            for gi, style in enumerate(edge_styles):
                edge_groups.append({
                    "class": style["class"],
                    "x": edge_flux_x_by_view[vi][gi],
                    "y": edge_flux_y_by_view[vi][gi],
                    "color": style["color"],
                    "width": style["width"],
                })

            views_data.append({
                "label": view.label,
                "rxn_nodes": rxn_node_list,
                "met_positions": [
                    {"x": flux_met_pos_by_view[vi][n][0],
                     "y": flux_met_pos_by_view[vi][n][1]}
                    for n in met_nodes
                ],
                "edge_groups": edge_groups,
            })

        # Metabolite node data
        met_node_list: list[dict] = []
        for n in met_nodes:
            consumers = [r for r, s in met_to_rxns.get(n, []) if s < 0]
            producers = [r for r, s in met_to_rxns.get(n, []) if s > 0]
            met_c = G.nodes[n].get("comp", "?")
            comp_cfg = compartments.get(met_c, {})
            met_node_list.append({
                "id": n,
                "name": G.nodes[n]["name"],
                "comp": met_c,
                "comp_label": comp_cfg.get("label", met_c),
                "comp_color": comp_cfg.get("color", "#8b949e"),
                "consumers": consumers[:8],
                "producers": producers[:8],
            })

        # Compartment hull vertices (data space — D3 renders paths via catmull-rom)
        compartment_list: list[dict] = []
        for comp_key, cfg_c in compartments.items():
            regular_here = [r for r in rxn_nodes
                            if G.nodes[r]["comp"] == comp_key
                            and G.nodes[r]["kind"] == "regular"
                            and r in pos]
            if len(regular_here) < 3:
                continue
            pts = np.array([[pos[r][0], pos[r][1]] for r in regular_here])
            try:
                hull = ConvexHull(pts)
                verts = pts[hull.vertices]
                centre = verts.mean(axis=0)
                verts = centre + (verts - centre) * 1.22
                compartment_list.append({
                    "key": comp_key,
                    "label": cfg_c["label"],
                    "color": cfg_c["color"],
                    "fill": cfg_c["fill"],
                    "hull_vertices": verts.tolist(),
                    "hull_tension": self._hull_tension,
                    "label_x": float(verts[:, 0].mean()),
                    "label_y": float(verts[:, 1].max()) + 0.2,
                })
            except Exception:
                pass

        # Data extent for D3 scales
        all_x = [pos[r][0] for r in rxn_nodes]
        all_y = [pos[r][1] for r in rxn_nodes]
        for view_pos in flux_met_pos_by_view:
            for n in met_nodes:
                all_x.append(view_pos[n][0])
                all_y.append(view_pos[n][1])
        pad = 1.0
        x_range = [min(all_x) - pad, max(all_x) + pad] if all_x else [-12.0, 8.0]
        y_range = [min(all_y) - pad, max(all_y) + pad] if all_y else [-9.0, 7.0]

        model_name = (mdl.id or "Model") if mdl else "Model"
        met_flux_x_by_view = [
            [flux_met_pos_by_view[vi][n][0] for n in met_nodes]
            for vi in range(n_views)
        ]
        met_flux_y_by_view = [
            [flux_met_pos_by_view[vi][n][1] for n in met_nodes]
            for vi in range(n_views)
        ]

        d3_data = {
            "meta": {
                "model_name": model_name,
                "x_range": x_range,
                "y_range": y_range,
                "height": 940,
                "abs_max_flux": float(abs_max_flux),
                "density_scale": float(ds),
            },
            "compartments": compartment_list,
            "view_labels": [v.label for v in views],
            "views": views_data,
            "met_nodes": met_node_list,
            "interactivity": {
                "view_labels": [v.label for v in views],
                "rxn_ids": rxn_nodes,
                "rxn_names": [mdl.reactions.get_by_id(r).name or r for r in rxn_nodes],
                "rxn_x": [pos[r][0] for r in rxn_nodes],
                "rxn_y": [pos[r][1] for r in rxn_nodes],
                "met_ids": met_nodes,
                "met_names": [G.nodes[n].get("name", n) for n in met_nodes],
                "met_flux_x": met_flux_x_by_view,
                "met_flux_y": met_flux_y_by_view,
            },
        }

        self._figure_json = json.dumps(d3_data, cls=_PlotlyJSONEncoder)

    # -- Frontend ESM --

    _esm = resources.files("cobramod.static").joinpath("flux_network.mjs").read_text()


class _PlotlyJSONEncoder(json.JSONEncoder):
    """Handle numpy types in JSON serialization."""
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)
