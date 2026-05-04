"""Layout and routing primitives for the FLoV widget.

Includes Fruchterman-Reingold spring layout, region scaling, hull-curve
sampling, A* grid routing, polyline simplification, and station branch
hub detection.

Split out of ``flux_network.py`` so the main file holds the widget class only.
"""
from __future__ import annotations

import heapq
from dataclasses import dataclass, field
from typing import Optional

import cobra
import networkx as nx
import numpy as np


# ==========================================================================
# DATACLASSES
# ==========================================================================

@dataclass
class HubLine:
    """One coloured line on a station route (≈ one metabolite class)."""
    klass: str
    rxn_ids: list[str]
    spine_offset_idx: int = 0


@dataclass
class StationPair:
    """A pair of hubs (one per compartment) connected by a routed spine."""
    pair_id: str
    comp_a: str
    comp_b: str
    lines: list[HubLine]
    anchor_a: tuple[float, float] = (0.0, 0.0)
    anchor_b: tuple[float, float] = (0.0, 0.0)
    tangent_a: tuple[float, float] = (0.0, 1.0)
    tangent_b: tuple[float, float] = (0.0, 1.0)
    spine: list[tuple[float, float]] = field(default_factory=list)


# ==========================================================================
# FRUCHTERMAN-REINGOLD SPRING LAYOUT
# ==========================================================================

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
    """
    nodes = list(G.nodes())
    n = len(nodes)
    if n == 0:
        return {}
    if n == 1:
        return {nodes[0]: (0.0, 0.0)}

    rng = np.random.default_rng(seed)
    pos = rng.random((n, 2)).astype(float) * 2.0 - 1.0

    if k is None:
        k = float(np.sqrt(1.0 / n))
    k2 = k * k

    node_idx = {v: i for i, v in enumerate(nodes)}
    if G.number_of_edges():
        u_arr = np.array([node_idx[u] for u, _ in G.edges()], dtype=np.intp)
        v_arr = np.array([node_idx[v] for _, v in G.edges()], dtype=np.intp)
    else:
        u_arr = v_arr = np.empty(0, dtype=np.intp)

    effective_iters = max(10, min(iterations, max(10, 10_000 // n)))

    t = max(0.1, 0.1 * np.sqrt(n))
    dt = t / (effective_iters + 1)

    delta  = np.empty((n, n, 2), dtype=float)
    dist2  = np.empty((n, n), dtype=float)
    scalar = np.empty((n, n), dtype=float)
    disp   = np.empty((n, 2),  dtype=float)

    for _ in range(effective_iters):
        np.subtract(pos[:, np.newaxis, :], pos[np.newaxis, :, :], out=delta)
        np.einsum("ijk,ijk->ij", delta, delta, out=dist2)
        np.fill_diagonal(dist2, 1.0)
        dist = np.sqrt(dist2)

        np.multiply(dist2, dist, out=scalar)
        np.divide(k2, scalar, out=scalar)
        np.fill_diagonal(scalar, 0.0)

        np.einsum("ij,ijk->ik", scalar, delta, out=disp)

        if u_arr.size:
            e_delta = pos[u_arr] - pos[v_arr]
            e_dist = np.hypot(e_delta[:, 0], e_delta[:, 1])
            e_dist = np.maximum(e_dist, 1e-10)
            att_mag = e_dist / k
            att_disp = e_delta / e_dist[:, np.newaxis] * att_mag[:, np.newaxis]
            disp[:, 0] -= np.bincount(u_arr, weights=att_disp[:, 0], minlength=n)
            disp[:, 0] += np.bincount(v_arr, weights=att_disp[:, 0], minlength=n)
            disp[:, 1] -= np.bincount(u_arr, weights=att_disp[:, 1], minlength=n)
            disp[:, 1] += np.bincount(v_arr, weights=att_disp[:, 1], minlength=n)

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


# ==========================================================================
# GRID + HULL SAMPLING HELPERS
# ==========================================================================

def _grid_snap(x: float, y: float, step: float) -> tuple[float, float]:
    return (round(x / step) * step, round(y / step) * step)


def _hull_curve_samples(
    verts: np.ndarray, tension: float = 0.4, n_per_seg: int = 12,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample (point, unit-tangent) along the smoothed Catmull-Rom hull.

    Mirrors the Catmull-Rom hull path rendered by the frontend so anchors land
    exactly on the rendered curve, not the raw polygon.
    """
    n = len(verts)
    if n < 3:
        return np.zeros((0, 2)), np.zeros((0, 2))
    pts: list[np.ndarray] = []
    tans: list[np.ndarray] = []
    ts = np.linspace(0.0, 1.0, n_per_seg, endpoint=False)
    for i in range(n):
        p0 = verts[(i - 1) % n]
        p1 = verts[i]
        p2 = verts[(i + 1) % n]
        p3 = verts[(i + 2) % n]
        cp1 = p1 + tension * (p2 - p0) / 3.0
        cp2 = p2 - tension * (p3 - p1) / 3.0
        for t in ts:
            mt = 1.0 - t
            pt = (mt ** 3) * p1 + 3 * (mt ** 2) * t * cp1 \
                + 3 * mt * (t ** 2) * cp2 + (t ** 3) * p2
            tan = (-3 * (mt ** 2) * p1
                   + 3 * (mt ** 2 - 2 * mt * t) * cp1
                   + 3 * (2 * mt * t - t ** 2) * cp2
                   + 3 * (t ** 2) * p2)
            pts.append(pt)
            tans.append(tan)
    pts_arr = np.asarray(pts)
    tans_arr = np.asarray(tans)
    norms = np.linalg.norm(tans_arr, axis=1, keepdims=True)
    norms[norms < 1e-9] = 1.0
    tans_arr = tans_arr / norms
    return pts_arr, tans_arr


def _project_to_curve(
    target: tuple[float, float], samples: np.ndarray,
) -> int:
    """Index of the curve sample closest to ``target``."""
    if len(samples) == 0:
        return -1
    diffs = samples - np.asarray(target)
    return int(np.argmin((diffs ** 2).sum(axis=1)))


def _polygon_mask(
    polygon: np.ndarray, gx: np.ndarray, gy: np.ndarray,
) -> np.ndarray:
    """Vectorised ray-casting point-in-polygon over a meshgrid.

    ``polygon`` is an Nx2 array of vertices (open ring is fine — last edge
    closes the polygon implicitly). ``gx``/``gy`` are 2-D meshgrid arrays.
    Returns a bool mask the same shape as ``gx``.
    """
    inside = np.zeros(gx.shape, dtype=bool)
    n = len(polygon)
    if n < 3:
        return inside
    j = n - 1
    eps = 1e-30
    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        cond1 = (yi > gy) != (yj > gy)
        cond2 = gx < (xj - xi) * (gy - yi) / (yj - yi + eps) + xi
        inside ^= (cond1 & cond2)
        j = i
    return inside


# ==========================================================================
# A* GRID ROUTING
# ==========================================================================

# 8 grid headings (drow, dcol) — N, NE, E, SE, S, SW, W, NW
_GRID_DIRS: tuple[tuple[int, int], ...] = (
    (-1, 0), (-1, 1), (0, 1), (1, 1),
    (1, 0), (1, -1), (0, -1), (-1, -1),
)
_DIR_COSTS: tuple[float, ...] = tuple(
    (dr * dr + dc * dc) ** 0.5 for dr, dc in _GRID_DIRS
)
_BEND_PENALTY = 0.4


def _a_star_grid(
    start: tuple[int, int],
    end: tuple[int, int],
    forbidden: np.ndarray,
    cell_penalty: Optional[np.ndarray] = None,
) -> Optional[list[tuple[int, int]]]:
    """8-direction A* with bend penalty. Returns cell list start→end or None.

    Heuristic: octile distance (admissible for 8-direction grids).
    Bend penalty is added on every direction change, which collapses
    stairway routes into single straights or one-bend Ls.
    """
    H, W = forbidden.shape
    if not (0 <= start[0] < H and 0 <= start[1] < W):
        return None
    if not (0 <= end[0] < H and 0 <= end[1] < W):
        return None
    if forbidden[start] or forbidden[end]:
        return None
    if cell_penalty is not None and cell_penalty.shape != forbidden.shape:
        raise ValueError("cell_penalty must match forbidden.shape")

    def heuristic(c: tuple[int, int]) -> float:
        dx = abs(c[1] - end[1]); dy = abs(c[0] - end[0])
        return (max(dx, dy) - min(dx, dy)) + (2 ** 0.5) * min(dx, dy)

    g_score: dict[tuple[tuple[int, int], int], float] = {(start, -1): 0.0}
    came_from: dict[tuple[tuple[int, int], int], tuple[tuple[int, int], int]] = {}
    heap: list[tuple[float, float, int, tuple[int, int], int]] = []
    counter = 0
    heapq.heappush(heap, (heuristic(start), 0.0, counter, start, -1))
    closed: set[tuple[tuple[int, int], int]] = set()

    while heap:
        _, g, _, cell, prev_di = heapq.heappop(heap)
        state = (cell, prev_di)
        if state in closed:
            continue
        closed.add(state)
        if cell == end:
            path = [cell]
            while state in came_from:
                state = came_from[state]
                path.append(state[0])
            return list(reversed(path))
        r, c = cell
        for di, (dr, dc) in enumerate(_GRID_DIRS):
            nr, nc = r + dr, c + dc
            if not (0 <= nr < H and 0 <= nc < W):
                continue
            if forbidden[nr, nc] and (nr, nc) != end:
                continue
            ng = g + _DIR_COSTS[di]
            if cell_penalty is not None:
                ng += float(cell_penalty[nr, nc])
            if prev_di != -1 and prev_di != di:
                ng += _BEND_PENALTY
            ns = ((nr, nc), di)
            if ng < g_score.get(ns, float("inf")):
                g_score[ns] = ng
                came_from[ns] = state
                counter += 1
                heapq.heappush(
                    heap, (ng + heuristic((nr, nc)), ng, counter, (nr, nc), di),
                )
    return None


# ==========================================================================
# POLYLINE HELPERS
# ==========================================================================

def _simplify_polyline(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """Drop collinear interior vertices — keeps only direction-change points."""
    if len(pts) <= 2:
        return list(pts)
    out = [pts[0]]
    for i in range(1, len(pts) - 1):
        ax, ay = pts[i - 1]
        bx, by = pts[i]
        cx, cy = pts[i + 1]
        if abs((bx - ax) * (cy - by) - (by - ay) * (cx - bx)) > 1e-9:
            out.append(pts[i])
    out.append(pts[-1])
    return out


# ==========================================================================
# BRANCH / MERGE HUB DETECTION
# ==========================================================================

def _detect_branch_hubs(
    stations: list[StationPair], step: float,
) -> list[dict]:
    """Find cells where two spines sharing an anchor diverge.

    Returns a list of ``{"x", "y", "branches": [pair_id, ...]}`` dicts.
    The first cell at which the spines part ways becomes a small branch
    pill in the rendered figure.
    """
    if not stations:
        return []
    by_anchor: dict[tuple[float, float], list[tuple[StationPair, list[tuple[float, float]]]]] = {}
    for st in stations:
        if not st.spine:
            continue
        snapped_a = _grid_snap(st.spine[0][0], st.spine[0][1], step)
        snapped_b = _grid_snap(st.spine[-1][0], st.spine[-1][1], step)
        by_anchor.setdefault(snapped_a, []).append((st, list(st.spine)))
        by_anchor.setdefault(snapped_b, []).append((st, list(reversed(st.spine))))

    branch_pts: list[dict] = []
    seen: set[tuple[float, float]] = set()
    for anchor, group in by_anchor.items():
        if len(group) < 2:
            continue
        max_len = min(len(s[1]) for s in group)
        diverged_at = None
        for k in range(1, max_len):
            cells = {_grid_snap(s[1][k][0], s[1][k][1], step) for s in group}
            if len(cells) > 1:
                diverged_at = k
                break
        if diverged_at is None:
            continue
        bx, by = _grid_snap(
            group[0][1][diverged_at - 1][0],
            group[0][1][diverged_at - 1][1],
            step,
        )
        if (bx, by) in seen or (bx, by) == anchor:
            continue
        seen.add((bx, by))
        branch_pts.append({
            "x": float(bx), "y": float(by),
            "branches": [s[0].pair_id for s in group],
        })
    return branch_pts


def _nearest_open_cell(
    target_xy: tuple[float, float],
    allowed_cells: np.ndarray,
    allowed_xy: np.ndarray,
    reserved: set[tuple[int, int]],
) -> Optional[tuple[int, int]]:
    """Nearest allowed cell to ``target_xy`` that is not reserved."""
    if allowed_cells.size == 0:
        return None
    if reserved:
        free_mask = np.array(
            [(int(rc[0]), int(rc[1])) not in reserved for rc in allowed_cells],
            dtype=bool,
        )
        if not free_mask.any():
            return None
        free_cells = allowed_cells[free_mask]
        free_xy = allowed_xy[free_mask]
    else:
        free_cells = allowed_cells
        free_xy = allowed_xy
    deltas = free_xy - np.asarray(target_xy, dtype=float)
    idx = int(np.argmin(np.einsum("ij,ij->i", deltas, deltas)))
    row, col = free_cells[idx]
    return int(row), int(col)
