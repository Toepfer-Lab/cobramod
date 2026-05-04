"""Builder mixin for the FLoV widget.

Holds the three heavy computation methods that turn a built graph + views
into the JSON figure descriptor:

* ``_compute_stations`` — bucket non-regular reactions into station pairs and
  route their spines through the exchange region with A*.
* ``_compute_regular_rails`` — route in-compartment reaction↔metabolite
  edges on a per-compartment grid, assigning lane slots.
* ``_rebuild_figure`` — assemble per-view trace data and serialise to JSON
  for the frontend.

Split out of ``flux_network.py`` so the main file holds the widget class
plumbing only.
"""
from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from ._flux_config import (
    ACTIVE_EDGE_EPS,
    MAX_LINES_PER_PAIR,
    MAX_RXNS_PER_LINE,
    _LINE_KLASS_COLORS,
    _LINE_KLASS_INACTIVE,
    _LINE_KLASSES,
    _RAIL_ROUTE_PALETTE,
)
from ._flux_helpers import (
    _edge_bucket,
    _emit_progress,
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
    _detect_branch_hubs,
    _detect_route_hubs,
    _grid_snap,
    _nearest_open_cell,
    _polygon_mask,
    _project_to_curve,
    _simplify_polyline,
)


class _BuilderMixin:
    """Methods that turn (graph, layout, views) into the JSON figure descriptor.

    Mixed into ``FLoV`` — relies on attributes set by the host class:
    ``_cobra_model``, ``_cfg``, ``_G``, ``_pos``, ``_views``, ``_rxn_nodes``,
    ``_met_nodes``, ``_met_to_rxns``, ``_rxn_kind_map``, ``_rxn_substrates``,
    ``_rxn_products``, ``_compartment_hulls``, ``_compartment_curve_samples``,
    ``_compartment_curve_tangents``, ``_effective_region``, plus the various
    user-facing flags like ``_drop_no_data`` / ``_stoich_flux`` / etc.
    """

    # -- Internal: station / railway builder --

    def _compute_stations(
        self, active_rxn_set: set[str], active_compartments: set[str],
    ) -> tuple[list[StationPair], list[dict], float, list[dict]]:
        """Build station list, run A* routing, detect branch hubs.

        Returns ``(stations, branch_hubs, grid_step, inner_edges)``.
        ``inner_edges`` are the rxn→met edges of non-regular reactions —
        kept as JSON so the frontend can toggle them on demand.
        """
        mdl = self._cobra_model
        cfg = self._cfg
        exch = cfg["exchange_key"]
        compartments = cfg["compartments"]
        G = self._G
        pos = self._pos
        if mdl is None or G is None or not pos:
            return [], [], 0.5, []

        # ── 1. Bucket reactions by (compA, compB) → klass → rxn_ids ──
        pair_buckets: dict[tuple[str, str], dict[str, list[str]]] = {}
        for rxn in mdl.reactions:
            if rxn.id not in active_rxn_set:
                continue
            kind = self._rxn_kind_map.get(rxn.id, "regular")
            if kind not in ("transport", "import", "export", "biomass"):
                continue
            non_exch = {_met_comp(m, cfg) for m in rxn.metabolites} - {exch}
            pk = _pair_key(kind, non_exch, cfg)
            if pk is None:
                continue
            if pk[0] not in self._compartment_hulls and pk[0] != exch:
                continue
            if pk[1] not in self._compartment_hulls and pk[1] != exch:
                continue
            klass = _met_class(rxn, cfg)
            pair_buckets.setdefault(pk, {}).setdefault(klass, []).append(rxn.id)

        if not pair_buckets:
            return [], [], 0.5, []

        # ── 2. Cap to MAX_LINES_PER_PAIR per pair, MAX_RXNS_PER_LINE per line ──
        stations: list[StationPair] = []
        for (a, b), klass_map in pair_buckets.items():
            ordered_lines: list[HubLine] = []
            for ki, klass in enumerate(_LINE_KLASSES):
                rxn_ids = klass_map.get(klass, [])
                if not rxn_ids:
                    continue
                ordered_lines.append(HubLine(
                    klass=klass,
                    rxn_ids=rxn_ids[:MAX_RXNS_PER_LINE]
                            if len(rxn_ids) > MAX_RXNS_PER_LINE else list(rxn_ids),
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

        # ── 3. Anchor each station on its compartment's smoothed hull curve ──
        all_x = [p[0] for p in pos.values()]
        all_y = [p[1] for p in pos.values()]
        x_min, x_max = (min(all_x), max(all_x)) if all_x else (-12.0, 8.0)
        y_min, y_max = (min(all_y), max(all_y)) if all_y else (-9.0, 7.0)
        pad = 2.0
        x_min -= pad; x_max += pad; y_min -= pad; y_max += pad
        span = max(x_max - x_min, y_max - y_min)
        step = max(0.4, span / 50.0)

        def _comp_centroid(ck: str) -> tuple[float, float]:
            if ck in self._compartment_hulls:
                v = self._compartment_hulls[ck]
                return float(v[:, 0].mean()), float(v[:, 1].mean())
            if ck in compartments:
                x0, x1, y0, y1 = compartments[ck]["region"]
                return (x0 + x1) / 2.0, (y0 + y1) / 2.0
            return (x_min + x_max) / 2.0, (y_min + y_max) / 2.0

        def _exch_anchor(other_centroid: tuple[float, float]) -> tuple[
            tuple[float, float], tuple[float, float],
        ]:
            """Synthesise an anchor for the exchange compartment along the
            edge of the global bounding box closest to *other_centroid*."""
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

        for st in stations:
            cb_x, cb_y = _comp_centroid(st.comp_b)
            ca_x, ca_y = _comp_centroid(st.comp_a)
            # Anchor A
            if st.comp_a in self._compartment_curve_samples:
                samples = self._compartment_curve_samples[st.comp_a]
                tans = self._compartment_curve_tangents[st.comp_a]
                idx = _project_to_curve((cb_x, cb_y), samples)
                ax, ay = samples[idx]
                tx, ty = tans[idx]
            else:
                (ax, ay), (tx, ty) = _exch_anchor((ca_x, ca_y))
            sx, sy = _grid_snap(ax, ay, step)
            st.anchor_a = (float(sx), float(sy))
            st.tangent_a = (float(tx), float(ty))
            # Anchor B
            if st.comp_b in self._compartment_curve_samples:
                samples = self._compartment_curve_samples[st.comp_b]
                tans = self._compartment_curve_tangents[st.comp_b]
                idx = _project_to_curve((ca_x, ca_y), samples)
                bx, by = samples[idx]
                tx2, ty2 = tans[idx]
            else:
                (bx, by), (tx2, ty2) = _exch_anchor((cb_x, cb_y))
            sx, sy = _grid_snap(bx, by, step)
            st.anchor_b = (float(sx), float(sy))
            st.tangent_b = (float(tx2), float(ty2))

        # ── 4. De-overlap multiple anchors at the same compartment ──
        for ck, samples in self._compartment_curve_samples.items():
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
                    tan = self._compartment_curve_tangents[ck][idx]
                    sx, sy = _grid_snap(pt[0], pt[1], step)
                    if side == "a":
                        st.anchor_a = (float(sx), float(sy))
                        st.tangent_a = (float(tan[0]), float(tan[1]))
                    else:
                        st.anchor_b = (float(sx), float(sy))
                        st.tangent_b = (float(tan[0]), float(tan[1]))
                prev_idx = idx

        # ── 5. Build forbidden bitmap from compartment hulls ──
        n_cols = int(np.ceil((x_max - x_min) / step)) + 1
        n_rows = int(np.ceil((y_max - y_min) / step)) + 1
        cols = np.arange(n_cols)
        rows = np.arange(n_rows)
        gx, gy = np.meshgrid(x_min + cols * step, y_min + rows * step, indexing="xy")
        forbidden = np.zeros((n_rows, n_cols), dtype=bool)
        for ck, verts in self._compartment_hulls.items():
            if self._drop_no_data and ck not in active_compartments:
                continue
            forbidden |= _polygon_mask(verts, gx, gy)

        def to_cell(p: tuple[float, float]) -> tuple[int, int]:
            r = int(round((p[1] - y_min) / step))
            c = int(round((p[0] - x_min) / step))
            r = max(0, min(n_rows - 1, r))
            c = max(0, min(n_cols - 1, c))
            return r, c

        def to_xy(cell: tuple[int, int]) -> tuple[float, float]:
            return (x_min + cell[1] * step, y_min + cell[0] * step)

        # ── 6. Route every station spine. ──
        for st in stations:
            sa = to_cell(st.anchor_a)
            sb = to_cell(st.anchor_b)
            saved = bool(forbidden[sa]), bool(forbidden[sb])
            forbidden[sa] = False
            forbidden[sb] = False
            cells = _a_star_grid(sa, sb, forbidden)
            forbidden[sa] = saved[0]
            forbidden[sb] = saved[1]
            if not cells:
                st.spine = [st.anchor_a, st.anchor_b]
                continue
            poly = [to_xy(c) for c in cells]
            poly[0] = st.anchor_a
            poly[-1] = st.anchor_b
            st.spine = _simplify_polyline(poly)

        # ── 7. Branch hubs ──
        branch_hubs = _detect_branch_hubs(stations, step)

        # ── 8. Inner edges ──
        inner_edges: list[dict] = []
        non_regular_active = {
            r for r in active_rxn_set
            if self._rxn_kind_map.get(r) in ("transport", "import", "export", "biomass")
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

        return stations, branch_hubs, step, inner_edges

    def _compute_regular_rails(
        self,
        rxn_nodes: list[str],
        met_nodes: list[str],
        views: list,
    ) -> dict:
        """Route regular in-compartment reaction↔metabolite links on a grid."""
        G = self._G
        pos = self._pos
        if G is None or not pos:
            return {
                "routes": [],
                "hubs": [],
                "guides": [],
                "rxn_pos": {r: pos.get(r, (0.0, 0.0)) for r in rxn_nodes},
                "met_pos": {m: pos.get(m, (0.0, 0.0)) for m in met_nodes},
                "routed_rxns": set(),
                "routed_mets": set(),
                "grid_step": 0.4,
            }

        edge_records: dict[str, list[dict]] = defaultdict(list)
        node_degree: dict[str, int] = defaultdict(int)
        flux_priority: dict[str, float] = defaultdict(float)
        for view in views:
            for rxn_id, flux in view.flux_dict.items():
                if _has_flux_data(float(flux)):
                    flux_priority[rxn_id] += abs(float(flux))

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
            if self._rxn_kind_map.get(rxn_id, "regular") != "regular":
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

        if not edge_records:
            return {
                "routes": [],
                "hubs": [],
                "guides": [],
                "rxn_pos": {r: pos.get(r, (0.0, 0.0)) for r in rxn_nodes},
                "met_pos": {m: pos.get(m, (0.0, 0.0)) for m in met_nodes},
                "routed_rxns": set(),
                "routed_mets": set(),
                "grid_step": 0.4,
            }

        routed_rxn_pos: dict[str, tuple[float, float]] = {}
        routed_met_pos: dict[str, tuple[float, float]] = {}
        routed_rxns: set[str] = set()
        routed_mets: set[str] = set()
        routes_geom: list[dict] = []
        route_hubs: list[dict] = []
        guide_links: list[dict] = []
        guide_seen: set[tuple[str, str]] = set()
        step_values: list[float] = []
        comp_track_xy: dict[str, np.ndarray] = {}

        for comp, records in edge_records.items():
            node_ids = sorted(
                {rec["rxn_id"] for rec in records} | {rec["met_id"] for rec in records},
                key=lambda nid: (-node_degree.get(nid, 0), nid),
            )
            if not node_ids:
                continue
            verts = self._compartment_hulls.get(comp)
            if verts is None or len(verts) < 3:
                x0, x1, y0, y1 = self._effective_region(comp)
                verts = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
            x_span = float(max(verts[:, 0].max() - verts[:, 0].min(), 1.0))
            y_span = float(max(verts[:, 1].max() - verts[:, 1].min(), 1.0))
            area = max(x_span * y_span, 1.0)
            step = float(np.clip((area / max(len(node_ids) * 7.0, 64.0)) ** 0.5, 0.16, 0.45))
            step_values.append(step)
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
            if int(allowed.sum()) < len(node_ids):
                allowed[:] = True

            allowed_cells = np.argwhere(allowed)
            allowed_xy = np.column_stack(
                (
                    x_min + allowed_cells[:, 1] * step,
                    y_min + allowed_cells[:, 0] * step,
                )
            )

            def to_xy(cell: tuple[int, int]) -> tuple[float, float]:
                return (
                    float(x_min + cell[1] * step),
                    float(y_min + cell[0] * step),
                )

            node_to_cell: dict[str, tuple[int, int]] = {}
            reserved: set[tuple[int, int]] = set()
            for node_id in node_ids:
                cell = _nearest_open_cell(pos[node_id], allowed_cells, allowed_xy, reserved)
                if cell is None:
                    continue
                node_to_cell[node_id] = cell
                reserved.add(cell)
                anchor_xy = to_xy(cell)
                if node_id in rxn_set:
                    routed_rxn_pos[node_id] = anchor_xy
                    routed_rxns.add(node_id)
                if node_id in met_set:
                    routed_met_pos[node_id] = anchor_xy
                    routed_mets.add(node_id)

            blocked = ~allowed
            track_usage = np.zeros((n_rows, n_cols), dtype=np.int16)
            route_cells: dict[str, list[tuple[int, int]]] = {}
            comp_used_cells: set[tuple[int, int]] = set()
            records.sort(
                key=lambda rec: (
                    -float(rec["priority"]),
                    -node_degree.get(rec["rxn_id"], 0),
                    -node_degree.get(rec["met_id"], 0),
                    str(rec["line_id"]),
                )
            )
            comp_start_idx = len(routes_geom)
            rec_by_lid: dict[str, dict] = {}
            ordered_lids: list[str] = []
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

            def lane_slot_candidates(
                max_lanes: int, n_routes: int,
            ) -> list[float]:
                """Route-local lane slots, ordered from visually centred out."""
                span = max(max_lanes, n_routes, 1)
                if max_lanes <= 1:
                    slots = [0.0]
                    for k in range(span):
                        slots.extend([-(k + 0.5), k + 0.5])
                    return slots
                if max_lanes % 2 == 0:
                    slots = []
                    for k in range(span):
                        slots.extend([-(k + 0.5), k + 0.5])
                    slots.append(0.0)
                    return slots
                slots = [0.0]
                for k in range(1, span + 1):
                    slots.extend([-float(k), float(k)])
                return slots

            used_edge_slots: dict[tuple, set[float]] = defaultdict(set)
            route_slots: dict[str, float] = {}
            for lid in ordered_lids:
                edge_sigs = route_edge_sigs[lid]
                max_lanes = max(
                    (len(cell_edge_lines[edge]) for edge, _ in edge_sigs),
                    default=1,
                )
                chosen_slot = 0.0
                candidates = lane_slot_candidates(max_lanes, len(ordered_lids))
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

            node_cells = set(node_to_cell.values())

            for lid in ordered_lids:
                rec = rec_by_lid[lid]
                cells = route_cells[lid]
                if len(cells) < 2:
                    path_pts = [to_xy(cells[0])] if cells else []
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
                    path_pts = [to_xy(cells[i]) for i in keep_idx]
                    route_slot = route_slots.get(lid, 0.0)
                    seg_slots = [route_slot] * (len(keep_idx) - 1)
                routes_geom.append({
                    "line_id": rec["line_id"],
                    "rxn_id": rec["rxn_id"],
                    "met_id": rec["met_id"],
                    "comp": comp,
                    "stoich": float(rec["stoich"]),
                    "path": [[float(p[0]), float(p[1])] for p in path_pts],
                    "seg_slots": seg_slots,
                })

            comp_track_xy[comp] = np.asarray(
                [to_xy(cell) for cell in sorted(comp_used_cells)],
                dtype=float,
            )
            excluded_cells = {
                cell for node_id, cell in node_to_cell.items()
                if node_id in routed_rxns or node_id in routed_mets
            }
            route_hubs.extend(_detect_route_hubs(route_cells, excluded_cells, to_xy))

            comp_routes = routes_geom[comp_start_idx:]
            for ci, r in enumerate(comp_routes):
                r["base_color"] = _RAIL_ROUTE_PALETTE[ci % len(_RAIL_ROUTE_PALETTE)]

        for rxn_id in rxn_nodes:
            routed_rxn_pos.setdefault(rxn_id, pos[rxn_id])
        for met_id in met_nodes:
            routed_met_pos.setdefault(met_id, pos[met_id])

        for node_id in rxn_nodes + met_nodes:
            comp = G.nodes[node_id].get("comp")
            track_xy = comp_track_xy.get(comp)
            if track_xy is None or track_xy.size == 0:
                continue
            routed = node_id in routed_rxns or node_id in routed_mets
            if routed:
                continue
            src = pos[node_id]
            deltas = track_xy - np.asarray(src, dtype=float)
            idx = int(np.argmin(np.einsum("ij,ij->i", deltas, deltas)))
            gx_, gy_ = track_xy[idx]
            node_type = "rxn" if node_id in rxn_set else "met"
            key = (node_type, node_id)
            if key in guide_seen:
                continue
            guide_seen.add(key)
            guide_links.append({
                "node_id": node_id,
                "node_type": node_type,
                "x": [float(src[0]), float(gx_)],
                "y": [float(src[1]), float(gy_)],
            })

        return {
            "routes": routes_geom,
            "hubs": route_hubs,
            "guides": guide_links,
            "rxn_pos": routed_rxn_pos,
            "met_pos": routed_met_pos,
            "routed_rxns": routed_rxns,
            "routed_mets": routed_mets,
            "grid_step": min(step_values) if step_values else 0.4,
        }

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

        # Drop-no-data filtering (on copies). Non-regular reactions render as
        # stations; regular reactions render as in-compartment rail nodes.
        rxn_nodes = [r for r in self._rxn_nodes
                     if self._rxn_kind_map.get(r, "regular") == "regular"]
        met_nodes = list(self._met_nodes)
        met_to_rxns = dict(self._met_to_rxns)

        # Active set across all views.
        active_rxn_set: set[str] = set()
        for view in views:
            for rxn_id, flux in view.flux_dict.items():
                if _has_flux_data(float(flux)):
                    active_rxn_set.add(rxn_id)
        if not active_rxn_set:
            active_rxn_set = {r for r in self._rxn_nodes
                              if self._rxn_kind_map.get(r) != "regular"}

        active_compartments: set[str] = set()
        for rxn_id in active_rxn_set:
            if rxn_id in G:
                comp = G.nodes[rxn_id].get("comp")
                if comp in compartments:
                    active_compartments.add(comp)

        # Non-regular reactions are rendered as stations between compartments.
        # Include both station endpoints so a compartment with only active
        # transport/exchange flux does not lose its hull.
        for rxn in mdl.reactions:
            if rxn.id not in active_rxn_set:
                continue
            kind = self._rxn_kind_map.get(rxn.id, "regular")
            if kind not in ("transport", "import", "export", "biomass"):
                continue
            non_exch = {_met_comp(m, cfg) for m in rxn.metabolites} - {exchange_key}
            pk = _pair_key(kind, non_exch, cfg)
            if pk is None:
                continue
            active_compartments.update(c for c in pk if c in compartments)

        if self._drop_no_data:
            rxn_nodes = [r for r in rxn_nodes if r in active_rxn_set]
            regular_remaining = set(rxn_nodes)
            active_station_mets: set[str] = set()
            for rxn_id in active_rxn_set:
                if self._rxn_kind_map.get(rxn_id, "regular") == "regular":
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
            met_nodes = [m for m in met_nodes if any(
                r in regular_remaining for r, _ in met_to_rxns.get(m, [])
            ) or m in active_station_mets]

        ds = self._resolve_density_scale(len(rxn_nodes))
        rail_data = self._compute_regular_rails(rxn_nodes, met_nodes, views)

        # ── Stations & routes ──
        stations, branch_hubs, grid_step, inner_edges = self._compute_stations(
            active_rxn_set,
            active_compartments,
        )
        pair_rxn_counts = [
            sum(len(line.rxn_ids) for line in st.lines) for st in stations
        ]
        median_count = float(np.median(pair_rxn_counts)) if pair_rxn_counts else 1.0
        median_count = max(median_count, 1.0)
        ROUTE_WMIN = max(1.5, 2.0 * ds)
        ROUTE_WMAX = max(ROUTE_WMIN + 1.5, 7.0 * ds)

        def _hull_short_dim(comp_key: str) -> float:
            v = self._compartment_hulls.get(comp_key)
            if v is not None and len(v) >= 3:
                return float(min(v[:, 0].max() - v[:, 0].min(),
                                 v[:, 1].max() - v[:, 1].min()))
            return 4.0

        # Global abs_max across all views (consistent colorscale)
        all_valid: list[float] = [
            float(f)
            for view in views
            for f in view.flux_dict.values()
            if _has_flux_data(float(f))
        ]
        abs_max_flux = max(float(np.nanmax(np.abs(all_valid))) if all_valid else 1.0, 1e-9)
        all_route_valid: list[float] = []
        for view in views:
            for route in rail_data["routes"]:
                f = view.flux_dict.get(route["rxn_id"], np.nan)
                if not _has_flux_data(float(f)):
                    continue
                factor = abs(float(route["stoich"])) if self._stoich_flux else 1.0
                all_route_valid.append(abs(float(f)) * factor)
        abs_max_route_flux = max(
            float(np.nanmax(all_route_valid)) if all_route_valid else 1.0,
            1e-9,
        )

        _exch_cfg = compartments.get(exchange_key, {"color": "#8b949e", "label": exchange_key})
        rail_wmin = max(0.9, 1.4 * ds)
        rail_wmax = max(rail_wmin + 0.7, 4.4 * ds)

        # Per-view data
        views_data: list[dict] = []
        rxn_render_pos = rail_data["rxn_pos"]
        met_render_pos = rail_data["met_pos"]
        routed_mets = rail_data["routed_mets"]
        for view in views:
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
                f, sd = float(fluxes[i]), float(stds[i])
                is_valid = bool(valid[i])
                fill_color = _flux_to_hex(f, abs_max_flux) if is_valid else "#6e7681"
                f_str = f"{f:+.4f} mmol/gDW/h" if is_valid else "--- (no data)"
                sd_str = f"{sd:.4f}" if not np.isnan(sd) else ""
                display_name = rxn_o.name or r
                subs = ", ".join(self._rxn_substrates.get(r, [])[:5]) or "---"
                prds = ", ".join(self._rxn_products.get(r, [])[:5]) or "---"
                rxn_node_list.append({
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

            # ── Per-view station line widths/colours ──
            station_views_data: list[dict] = []
            for st in stations:
                line_view: list[dict] = []
                for ln in st.lines:
                    flux_sum = 0.0
                    for rid in ln.rxn_ids:
                        f = view.flux_dict.get(rid, np.nan)
                        if _has_flux_data(float(f)):
                            flux_sum += abs(float(f))
                    norm = float(np.clip(flux_sum / abs_max_flux, 0.0, 1.0))
                    width = ROUTE_WMIN + (ROUTE_WMAX - ROUTE_WMIN) * norm
                    color = (_LINE_KLASS_COLORS.get(ln.klass, "#ec4899")
                             if flux_sum > ACTIVE_EDGE_EPS
                             else _LINE_KLASS_INACTIVE)
                    line_view.append({
                        "klass": ln.klass,
                        "flux_sum": flux_sum,
                        "width": float(width),
                        "color": color,
                    })
                station_views_data.append({
                    "pair_id": st.pair_id,
                    "lines": line_view,
                })

            rail_views_data: list[dict] = []
            for route in rail_data["routes"]:
                f = view.flux_dict.get(route["rxn_id"], np.nan)
                if _has_flux_data(float(f)):
                    eff_flux = float(f) * (
                        abs(float(route["stoich"])) if self._stoich_flux else 1.0
                    )
                    _, mag_norm, sign = _edge_bucket(eff_flux, abs_max_route_flux)
                    width = rail_wmax
                    color = _modulate_route_color(
                        route.get("base_color", "#58a6ff"),
                        mag_norm,
                        True,
                    )
                else:
                    width = rail_wmax
                    color = "rgba(110,118,129,0.18)"
                rail_views_data.append({
                    "line_id": route["line_id"],
                    "width": float(width),
                    "color": color,
                })

            views_data.append({
                "label": view.label,
                "rxn_nodes": rxn_node_list,
                "rail_routes": rail_views_data,
                "stations": station_views_data,
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
                "x": met_render_pos[n][0],
                "y": met_render_pos[n][1],
                "comp": met_c,
                "comp_label": comp_cfg.get("label", met_c),
                "comp_color": comp_cfg.get("color", "#8b949e"),
                "display_kind": "stop" if n in routed_mets else "detached",
                "consumers": consumers[:8],
                "producers": producers[:8],
            })

        # Compartment hull vertices
        compartment_list: list[dict] = []
        for comp_key, cfg_c in compartments.items():
            if self._drop_no_data and comp_key not in active_compartments:
                continue
            verts = self._compartment_hulls.get(comp_key)
            if verts is None or len(verts) < 3:
                continue
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

        model_name = (mdl.id or "Model") if mdl else "Model"

        # ── Station geometry ──
        stations_geometry: list[dict] = []
        for st in stations:
            n_rxns_total = sum(len(line.rxn_ids) for line in st.lines)
            scale = float(np.log1p(n_rxns_total) / max(np.log1p(median_count), 1e-6))
            scale = float(np.clip(scale, 0.85, 1.4))
            local_size = min(_hull_short_dim(st.comp_a), _hull_short_dim(st.comp_b))
            pill_h = max(0.25, 0.10 * local_size * scale)
            pill_w = 0.35 * pill_h
            n_lines = len(st.lines)
            comp_a_cfg = compartments.get(st.comp_a, _exch_cfg)
            comp_b_cfg = compartments.get(st.comp_b, _exch_cfg)
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
                        "substrates": self._rxn_substrates.get(rid, [])[:5],
                        "products": self._rxn_products.get(rid, [])[:5],
                    })
                line_geom.append({
                    "klass": ln.klass,
                    "rxn_ids": ln.rxn_ids,
                    "rxn_summaries": rxn_summaries,
                    "spine_offset_idx": ln.spine_offset_idx,
                    "path": [[float(p[0]), float(p[1])] for p in st.spine],
                    "seg_slots": [lane_slot] * max(len(st.spine) - 1, 0),
                })
            stations_geometry.append({
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

        # Data extent for D3 scales
        all_x = [rxn_render_pos[r][0] for r in rxn_nodes] + [met_render_pos[n][0] for n in met_nodes]
        all_y = [rxn_render_pos[r][1] for r in rxn_nodes] + [met_render_pos[n][1] for n in met_nodes]
        for route in rail_data["routes"]:
            for px, py in route["path"]:
                all_x.append(float(px))
                all_y.append(float(py))
        for hub in rail_data["hubs"] + branch_hubs:
            all_x.append(float(hub["x"]))
            all_y.append(float(hub["y"]))
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

        d3_data = {
            "meta": {
                "model_name": model_name,
                "x_range": x_range,
                "y_range": y_range,
                "height": 940,
                "abs_max_flux": float(abs_max_flux),
                "density_scale": float(ds),
                "grid_step": float(min(grid_step, rail_data["grid_step"])),
            },
            "compartments": compartment_list,
            "view_labels": [v.label for v in views],
            "views": views_data,
            "met_nodes": met_node_list,
            "rail_routes": rail_data["routes"],
            "rail_hubs": rail_data["hubs"],
            "guide_links": rail_data["guides"],
            "stations": stations_geometry,
            "branch_hubs": branch_hubs,
            "inner_edges": inner_edges,
            "interactivity": {
                "view_labels": [v.label for v in views],
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

        self._figure_json = json.dumps(d3_data, cls=_PlotlyJSONEncoder)


# ==========================================================================
# JSON ENCODER
# ==========================================================================

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
