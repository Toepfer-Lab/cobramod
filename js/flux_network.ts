// D3-based renderer for the FLoV metabolic flux network anywidget
import * as d3 from 'd3';

// ── Interfaces ────────────────────────────────────────────────────────────────

interface HoverData {
  display_name: string; id: string; kind_badge: string;
  comp_label: string; pipeline_tag: string;
  flux_str: string; std_str: string;
  substrates: string; products: string; extra: string;
}
interface RxnNode {
  id: string; name: string; x: number; y: number;
  comp: string; kind: string; flux: number | null;
  fill_color: string; border_color: string;
  border_width: number; opacity: number; r: number;
  symbol: string; hover: HoverData;
}
interface RailRouteView {
  line_id: string; width: number; color: string;
}
interface StationLineView {
  klass: string;
  flux_sum: number;
  width: number;
  color: string;
}
interface StationViewEntry {
  pair_id: string;
  lines: StationLineView[];
}
interface ViewData {
  label: string; rxn_nodes: RxnNode[];
  rail_routes: RailRouteView[];
  stations: StationViewEntry[];
}
interface RxnSummary {
  id: string; name: string; flux_per_view: string[];
}
interface RailRouteGeom {
  line_id: string; rxn_id: string; met_id: string;
  comp: string; stoich: number;
  base_color: string;
  seg_slots: number[];
  path: [number, number][];
}
interface StationLineGeom {
  klass: string;
  rxn_ids: string[];
  rxn_summaries: RxnSummary[];
  spine_offset_idx: number;
  path: [number, number][];
}
interface StationGeom {
  pair_id: string;
  comp_a: string; comp_b: string;
  comp_a_label: string; comp_b_label: string;
  comp_a_color: string; comp_b_color: string;
  anchor_a: [number, number]; anchor_b: [number, number];
  tangent_a: [number, number]; tangent_b: [number, number];
  pill_height: number; pill_width: number;
  n_rxns: number;
  lines: StationLineGeom[];
}
interface BranchHub {
  x: number; y: number; branches: string[];
}
interface InnerEdge {
  rxn_id: string; met_id: string;
  x: [number, number]; y: [number, number];
}
interface RailHub {
  x: number; y: number; role: string; routes: string[];
}
interface GuideLink {
  node_id: string; node_type: 'rxn' | 'met';
  x: [number, number]; y: [number, number];
}
interface MetNode {
  id: string; name: string; x: number; y: number;
  comp: string; comp_label: string; comp_color: string;
  display_kind: 'stop' | 'detached';
  consumers: string[]; producers: string[];
}
interface Compartment {
  key: string; label: string; color: string; fill: string;
  hull_vertices: [number, number][];
  hull_tension: number; label_x: number; label_y: number;
}
interface Meta {
  model_name: string; x_range: [number, number]; y_range: [number, number];
  height: number; abs_max_flux: number; density_scale: number;
  grid_step: number;
}
interface Interactivity {
  view_labels: string[];
  rxn_ids: string[]; rxn_names: string[]; rxn_x: number[]; rxn_y: number[];
  met_ids: string[]; met_names: string[]; met_x: number[]; met_y: number[];
}
interface D3FigureData {
  meta: Meta; compartments: Compartment[];
  view_labels: string[]; views: ViewData[];
  met_nodes: MetNode[]; interactivity: Interactivity;
  rail_routes: RailRouteGeom[]; rail_hubs: RailHub[]; guide_links: GuideLink[];
  stations: StationGeom[];
  branch_hubs: BranchHub[];
  inner_edges: InnerEdge[];
}

// ── CSS ───────────────────────────────────────────────────────────────────────

const CSS = `
.flov-wrapper{background:#0d1117;position:relative;font-family:'Inter',Arial,sans-serif;overflow:hidden;}
.flov-toolbar{display:flex;flex-wrap:wrap;gap:6px;padding:8px 12px;background:#161b22;border-bottom:1px solid #30363d;align-items:flex-end;}
.flov-btn-group{display:flex;gap:3px;align-items:flex-end;}
.flov-btn-label{font:10px/1 'Inter',Arial,sans-serif;color:#7d8590;padding:0 4px 5px 0;}
.flov-btn{background:#21262d;color:#e6edf3;border:1px solid #30363d;border-radius:6px;padding:5px 13px;font:11px 'Inter',Arial,sans-serif;cursor:pointer;transition:background 80ms,border-color 80ms;white-space:nowrap;}
.flov-btn:hover{background:#30363d;border-color:#8b949e;}
.flov-btn.active{background:#1f6feb;border-color:#1f6feb;color:#fff;}
.flov-svg{display:block;}
.flov-tt{position:absolute;pointer-events:none;z-index:1000;background:#1c2128;border-radius:8px;padding:10px 14px;font:12px/1.65 'Inter',Arial,sans-serif;color:#e6edf3;max-width:290px;box-shadow:0 6px 24px rgba(0,0,0,0.75);display:none;border:1.5px solid #30363d;}
.flov-tt b{color:#fff;}
.flov-tt .ft-flux{font-family:monospace;color:#ffa657;}
.flov-tt .ft-div{color:#30363d;}
.flov-tt .ft-comp{color:#8b949e;font-style:italic;}
.flov-tt .ft-badge{background:#21262d;border-radius:3px;padding:1px 5px;font-size:10px;color:#8b949e;letter-spacing:.04em;}
.flov-search-wrap{position:absolute;top:10px;right:10px;z-index:100;background:#161b22;border:1px solid #30363d;border-radius:6px;padding:5px 10px;width:230px;box-shadow:0 4px 14px rgba(0,0,0,0.55);}
.flov-search-input{background:transparent;border:none;outline:none;color:#e6edf3;width:100%;font:11px 'Inter',Arial,sans-serif;}
.flov-search-input::placeholder{color:#7d8590;}
.flov-search-results{display:none;max-height:184px;overflow-y:auto;margin-top:4px;}
.flov-search-result{padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #21262d;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;color:#e6edf3;}
.flov-search-result:hover{background:#21262d;}
@keyframes flov-pulse{0%,100%{stroke-opacity:1;stroke-width:3}50%{stroke-opacity:.35;stroke-width:7}}
.flov-hl{animation:flov-pulse 1.3s ease-in-out infinite;}
.flov-wrapper.light-mode{background:#f6f8fa;color:#24292f;}
.flov-wrapper.light-mode .flov-toolbar{background:#ffffff;border-bottom-color:#d0d7de;}
.flov-wrapper.light-mode .flov-btn-group{color:#24292f;}
.flov-wrapper.light-mode .flov-btn-label{color:#57606a;}
.flov-wrapper.light-mode .flov-btn{background:#f6f8fa;color:#24292f;border-color:#d0d7de;}
.flov-wrapper.light-mode .flov-btn:hover{background:#eef2f7;border-color:#8c959f;}
.flov-wrapper.light-mode .flov-btn.active{background:#0969da;border-color:#0969da;color:#ffffff;}
.flov-wrapper.light-mode .flov-tt{background:#ffffff;color:#24292f;border-color:#d0d7de;box-shadow:0 10px 28px rgba(31,35,40,0.16);}
.flov-wrapper.light-mode .flov-tt b{color:#24292f;}
.flov-wrapper.light-mode .flov-tt .ft-flux{color:#b35900;}
.flov-wrapper.light-mode .flov-tt .ft-div,.flov-wrapper.light-mode .flov-tt .ft-comp,.flov-wrapper.light-mode .flov-tt .ft-badge{color:#57606a;}
.flov-wrapper.light-mode .flov-tt .ft-badge{background:#f3f4f6;}
.flov-wrapper.light-mode .flov-search-wrap{background:#ffffff;border-color:#d0d7de;box-shadow:0 8px 24px rgba(31,35,40,0.12);}
.flov-wrapper.light-mode .flov-search-input{color:#24292f;}
.flov-wrapper.light-mode .flov-search-input::placeholder{color:#57606a;}
.flov-wrapper.light-mode .flov-search-result{color:#24292f;border-top-color:#d8dee4;}
.flov-wrapper.light-mode .flov-search-result:hover{background:#eef2f7;}
`;

let _cssInjected = false;
function injectCSS(): void {
  if (_cssInjected) return;
  const s = document.createElement('style');
  s.textContent = CSS;
  document.head.appendChild(s);
  _cssInjected = true;
}

// ── Symbol helpers ────────────────────────────────────────────────────────────

// Custom flat-top hexagon symbol (transport reactions)
const symbolHexagon: d3.SymbolType = {
  draw(ctx: { moveTo(x:number,y:number):void; lineTo(x:number,y:number):void; closePath():void }, size: number): void {
    const r = Math.sqrt(size / (3 * Math.sqrt(3) / 2));
    for (let i = 0; i < 6; i++) {
      const a = (i * Math.PI) / 3;
      i === 0 ? ctx.moveTo(r * Math.cos(a), r * Math.sin(a))
              : ctx.lineTo(r * Math.cos(a), r * Math.sin(a));
    }
    ctx.closePath();
  },
};

// r is the Plotly-style diameter (px); D3 symbol size = area in px²
function symArea(r: number): number { return Math.max(Math.PI * (r / 2) * (r / 2), 6); }

function symPath(symbol: string, r: number): string {
  const area = symArea(r);
  switch (symbol) {
    case 'triangle-up':   return d3.symbol(d3.symbolTriangle, area)(null) ?? '';
    case 'triangle-down': return d3.symbol(d3.symbolTriangle, area)(null) ?? '';
    case 'hexagon':       return d3.symbol(symbolHexagon,      area)(null) ?? '';
    case 'cross':         return d3.symbol(d3.symbolCross,    area * 0.7)(null) ?? '';
    default:              return d3.symbol(d3.symbolCircle,   area)(null) ?? '';
  }
}

function metPath(m: MetNode): string {
  if (m.display_kind === 'stop') {
    return d3.symbol(d3.symbolCircle, Math.PI * 20)(null) ?? '';
  }
  return d3.symbol(d3.symbolDiamond, Math.PI * 16)(null) ?? '';
}

// ── Path utilities ────────────────────────────────────────────────────────────

function hullPath(
  verts: [number, number][],
  xS: d3.ScaleLinear<number, number>, yS: d3.ScaleLinear<number, number>
): string {
  const gen = d3.line<[number, number]>()
    .x(d => xS(d[0])).y(d => yS(d[1]))
    .curve(d3.curveCatmullRomClosed.alpha(0.5));
  return gen(verts) ?? '';
}

// Polyline path through orthogonal/45° spine vertices (no curve smoothing —
// the bend penalty in the A* router has already collapsed stairways into
// single straights or one-bend Ls).
function polylinePath(
  pts: [number, number][],
  xS: d3.ScaleLinear<number, number>, yS: d3.ScaleLinear<number, number>
): string {
  if (!pts.length) return '';
  let d = `M${xS(pts[0][0]).toFixed(1)} ${yS(pts[0][1]).toFixed(1)}`;
  for (let i = 1; i < pts.length; i++) {
    d += ` L${xS(pts[i][0]).toFixed(1)} ${yS(pts[i][1]).toFixed(1)}`;
  }
  return d;
}

// Rail path with pixel-constant lane offset.
// Each segment is shifted right-perpendicularly by seg_slots[k] * LANE_W_PX
// divided by the current zoom scale k_zoom, so the screen-pixel gap stays
// constant as the user zooms. Miter intersections keep every segment truly
// parallel; sharp corners fall back to a bevel midpoint.
const LANE_W_PX = 3.5;
function railOffsetPath(
  pts: [number, number][],
  seg_slots: number[],
  xS: d3.ScaleLinear<number, number>,
  yS: d3.ScaleLinear<number, number>,
  k_zoom: number
): string {
  const n = pts.length;
  if (n < 2) return '';
  // SVG-pixel coordinates of the center-line path.
  const sx = pts.map(p => xS(p[0]));
  const sy = pts.map(p => yS(p[1]));
  // Per-segment unit directions and right-hand normals in SVG-pixel space.
  const ux: number[] = [], uy: number[] = [];
  const nx: number[] = [], ny: number[] = [];
  for (let k = 0; k < n - 1; k++) {
    const dx = sx[k + 1] - sx[k], dy = sy[k + 1] - sy[k];
    const len = Math.hypot(dx, dy) || 1;
    ux.push(dx / len); uy.push(dy / len);
    nx.push(dy / len); ny.push(-dx / len);  // right-hand ⊥
  }
  // Pixel offset per segment (constant in screen space regardless of zoom).
  const off = seg_slots.map(s => (s * LANE_W_PX) / k_zoom);
  // Build offset path vertices using miter joins at interior vertices.
  const MITER_CAP = (4 * LANE_W_PX) / k_zoom;
  const pts_out: [number, number][] = [];
  // First endpoint: offset perpendicular to first segment.
  pts_out.push([sx[0] + off[0] * nx[0], sy[0] + off[0] * ny[0]]);
  for (let k = 1; k < n - 1; k++) {
    // Offset base points on each adjacent segment at vertex k.
    const Ax = sx[k] + off[k-1] * nx[k-1], Ay = sy[k] + off[k-1] * ny[k-1];
    const Bx = sx[k] + off[k]   * nx[k],   By = sy[k] + off[k]   * ny[k];
    // Miter: intersect the two offset lines.
    const cross = ux[k-1] * uy[k] - uy[k-1] * ux[k];
    if (Math.abs(cross) < 1e-9) {
      // Collinear segments. If the lane offset changes (e.g. a sibling route
      // just terminated at this junction), emit two distinct corner points
      // so the line steps perpendicularly between lanes — a sharp jog at the
      // node instead of a diagonal across the whole segment.
      if (Math.abs(off[k] - off[k-1]) < 1e-6) {
        pts_out.push([(Ax + Bx) / 2, (Ay + By) / 2]);
      } else {
        pts_out.push([Ax, Ay]);
        pts_out.push([Bx, By]);
      }
      continue;
    } else {
      const ddx = Bx - Ax, ddy = By - Ay;
      const t = (ddx * uy[k] - ddy * ux[k]) / cross;
      const mx = Ax + t * ux[k-1], my = Ay + t * uy[k-1];
      const dist = Math.hypot(mx - sx[k], my - sy[k]);
      if (dist > MITER_CAP) {
        pts_out.push([(Ax + Bx) / 2, (Ay + By) / 2]);
      } else {
        pts_out.push([mx, my]);
      }
    }
  }
  // Last endpoint: offset perpendicular to last segment.
  pts_out.push([sx[n-1] + off[n-2] * nx[n-2], sy[n-1] + off[n-2] * ny[n-2]]);
  let d = `M${pts_out[0][0].toFixed(1)} ${pts_out[0][1].toFixed(1)}`;
  for (let i = 1; i < pts_out.length; i++) {
    d += ` L${pts_out[i][0].toFixed(1)} ${pts_out[i][1].toFixed(1)}`;
  }
  return d;
}

// ── Tooltip HTML ──────────────────────────────────────────────────────────────

function rxnTooltipHTML(h: HoverData): string {
  const badge = h.kind_badge ? `<span class="ft-badge">${h.kind_badge}</span>` : '';
  const id    = h.id ? ` <span style="color:#8b949e;font-size:11px">(${h.id})</span>` : '';
  const pipe  = h.pipeline_tag ? ` · ${h.pipeline_tag}` : '';
  const std   = h.std_str ? `<br><span style="color:#7d8590;font-size:10px">σ = ${h.std_str}</span>` : '';
  const extra = h.extra ? `<br>${h.extra}` : '';
  return (
    `<b>${h.display_name}</b>${id} ${badge}<br>` +
    `<span class="ft-div">──────────────────────</span><br>` +
    `<span class="ft-comp">${h.comp_label}${pipe}</span><br>` +
    `Flux: <span class="ft-flux">${h.flux_str}</span>${std}<br>` +
    `<span style="color:#8b949e">Sub:</span> ${h.substrates}<br>` +
    `<span style="color:#8b949e">Prd:</span> ${h.products}${extra}`
  );
}

function metTooltipHTML(m: MetNode): string {
  return (
    `<b>${m.name}</b><br>` +
    `<span class="ft-div">──────────────────────</span><br>` +
    `<span class="ft-comp">${m.comp_label}</span><br>` +
    `<span style="color:#8b949e">ID:</span> ${m.id}<br>` +
    `<span style="color:#8b949e">Consumed:</span> ${m.consumers.join(', ') || '—'}<br>` +
    `<span style="color:#8b949e">Produced:</span> ${m.producers.join(', ') || '—'}`
  );
}

// Per-class colour palette — must mirror _LINE_KLASS_COLORS in
// flux_network.py.  Used for the swatch shown next to each class header in
// the station tooltip so the on-route colour is immediately legible.
const KLASS_COLORS: Record<string, string> = {
  amino_acid: '#22c55e',
  sugar:      '#f59e0b',
  cofactor:   '#a855f7',
  inorganic:  '#06b6d4',
  other:      '#ec4899',
};

// Station tooltip: lists every reaction grouped under the station,
// sub-grouped by metabolite class with per-view fluxes.
function stationTooltipHTML(st: StationGeom): string {
  const head =
    `<b>${st.comp_a_label} ↔ ${st.comp_b_label}</b> ` +
    `<span class="ft-badge">${st.n_rxns} rxn${st.n_rxns === 1 ? '' : 's'}</span><br>` +
    `<span class="ft-div">──────────────────────</span><br>`;
  const lineRows = st.lines.map(ln => {
    const swatchColor = KLASS_COLORS[ln.klass] ?? '#6e7681';
    const swatch = `<span style="display:inline-block;width:9px;height:9px;border-radius:2px;background:${swatchColor};margin-right:5px;vertical-align:-1px"></span>`;
    const rxns = ln.rxn_summaries.slice(0, 5).map(r =>
      `<div style="margin-left:8px"><span style="color:#8b949e">${r.id}</span> ` +
      `${r.name && r.name !== r.id ? '· ' + r.name : ''}<br>` +
      `<span style="margin-left:10px;color:#7d8590;font-size:10px">${r.flux_per_view.join(' · ')}</span></div>`
    ).join('');
    const more = ln.rxn_summaries.length > 5
      ? `<div style="margin-left:8px;color:#7d8590;font-size:10px">…+${ln.rxn_summaries.length - 5} more</div>`
      : '';
    return `<div>${swatch}<b style="color:${swatchColor}">${ln.klass}</b></div>${rxns}${more}`;
  }).join('');
  return head + lineRows;
}

function branchTooltipHTML(b: BranchHub): string {
  return (
    `<b>Branch hub</b><br>` +
    `<span class="ft-div">──────────────────────</span><br>` +
    `<span style="color:#8b949e">Diverging routes:</span><br>` +
    b.branches.map(p => `<div style="margin-left:8px">${p.replace('|', ' ↔ ')}</div>`).join('')
  );
}

function railHubTooltipHTML(h: RailHub): string {
  return (
    `<b>${h.role}</b><br>` +
    `<span class="ft-div">──────────────────────</span><br>` +
    `<span style="color:#8b949e">Connected rails:</span><br>` +
    h.routes.map(r => `<div style="margin-left:8px">${r.replace('->', ' → ')}</div>`).join('')
  );
}

// ── Main render function ──────────────────────────────────────────────────────

function render({ model: S, el: I }: { model: any; el: HTMLElement }): void {
  injectCSS();

  const wrapper = document.createElement('div');
  wrapper.className = 'flov-wrapper';
  I.appendChild(wrapper);

  // Mutable state
  let fig: D3FigureData | null = null;
  let currentView = 0;
  let xS: d3.ScaleLinear<number, number> = d3.scaleLinear();
  let yS: d3.ScaleLinear<number, number> = d3.scaleLinear();
  let svgEl: SVGSVGElement | null = null;
  let zoomB: d3.ZoomBehavior<SVGSVGElement, unknown> | null = null;
  let viewBtns: HTMLButtonElement[] = [];
  let tooltipEl: HTMLDivElement | null = null;
  let resizeTimer: ReturnType<typeof setTimeout> | null = null;
  let gMetSel: d3.Selection<SVGGElement, unknown, null, undefined> | null = null;
  let gMetHitSel: d3.Selection<SVGPathElement, MetNode, SVGGElement, unknown> | null = null;
  let gRxnSel: d3.Selection<SVGPathElement, RxnNode, SVGGElement, unknown> | null = null;
  let gRxnHitSel: d3.Selection<SVGPathElement, RxnNode, SVGGElement, unknown> | null = null;
  let currentTransform: d3.ZoomTransform = d3.zoomIdentity;
  let rafId: number | null = null;
  let cachedWrapperRect: DOMRect | null = null;
  let lightMode = true;

  const MARGIN = { t: 10, r: 90, b: 10, l: 20 };
  const MET_HOVER_STROKE = 12;
  const RXN_HOVER_STROKE = 14;

  function metTransform(node: MetNode, transform: d3.ZoomTransform): string {
    const tx = xS(node.x).toFixed(1);
    const ty = yS(node.y).toFixed(1);
    const inv = (1 / transform.k).toFixed(6);
    return `matrix(${inv},0,0,${inv},${tx},${ty})`;
  }

  function updateMetZoom(transform: d3.ZoomTransform): void {
    if (!gMetSel) return;
    gMetSel.selectAll<SVGPathElement, MetNode>('path.met-shape')
      .attr('transform', d => metTransform(d, transform));
    gMetHitSel?.attr('transform', d => metTransform(d, transform));
  }

  function updateRailZoom(transform: d3.ZoomTransform): void {
    const zl = wrapper.querySelector<SVGGElement>('.zoom-layer');
    if (!zl || !fig) return;
    d3.select(zl).select('.g-rail')
      .selectAll<SVGPathElement, RailRouteGeom>('path.rail-line')
      .attr('d', d => railOffsetPath(d.path, d.seg_slots, xS, yS, transform.k));
  }

  function updateRxnZoom(transform: d3.ZoomTransform): void {
    if (!gRxnSel || !gRxnHitSel) return;
    gRxnSel.attr('transform', (d: RxnNode) => rxnTransform(d, transform));
    gRxnHitSel.attr('transform', (d: RxnNode) => rxnTransform(d, transform));
  }

  // ── Build full UI ──

  function buildUI(data: D3FigureData): void {
    fig = data;
    currentView = 0;
    currentTransform = d3.zoomIdentity;
    gMetSel = null;
    gMetHitSel = null;
    gRxnSel = null;
    gRxnHitSel = null;
    cachedWrapperRect = null;
    if (rafId !== null) { cancelAnimationFrame(rafId); rafId = null; }
    wrapper.innerHTML = '';
    tooltipEl = null;
    viewBtns = [];

    const svgW = Math.max(wrapper.clientWidth || 900, 600);
    const svgH = data.meta.height;

    xS = d3.scaleLinear().domain(data.meta.x_range).range([MARGIN.l, svgW - MARGIN.r]);
    yS = d3.scaleLinear().domain(data.meta.y_range).range([svgH - MARGIN.b, MARGIN.t]);

    wrapper.classList.toggle('light-mode', lightMode);
    buildToolbar(data, wrapper);
    const { svg, zoomLayer } = buildSVG(svgW, svgH);
    buildColorbar(wrapper, data.meta.abs_max_flux, svgH);
    buildSearch(wrapper, data, svg.node()!, svgW, svgH);
    buildTooltip(wrapper);
    drawAll(zoomLayer, data, svgW);

    // Zoom
    zoomB = d3.zoom<SVGSVGElement, unknown>()
      .scaleExtent([0.12, 40])
      .on('zoom', (e: d3.D3ZoomEvent<SVGSVGElement, unknown>) => {
        currentTransform = e.transform;
        zoomLayer.attr('transform', e.transform.toString());
        if (rafId !== null) cancelAnimationFrame(rafId);
        rafId = requestAnimationFrame(() => {
          rafId = null;
          updateRxnZoom(currentTransform);
          updateMetZoom(currentTransform);
          updateRailZoom(currentTransform);
        });
      });
    svg.call(zoomB);
    svg.on('dblclick.zoom', () => {
      svg.transition().duration(380).call(zoomB!.transform, d3.zoomIdentity);
    });
  }

  // ── SVG skeleton ──

  function buildSVG(
    svgW: number, svgH: number
  ): { svg: d3.Selection<SVGSVGElement, unknown, null, undefined>; zoomLayer: d3.Selection<SVGGElement, unknown, null, undefined> } {
    const svg = d3.select(wrapper)
      .append<SVGSVGElement>('svg')
      .attr('class', 'flov-svg')
      .attr('width', svgW)
      .attr('height', svgH)
      .style('background', lightMode ? '#f6f8fa' : '#0d1117');
    svgEl = svg.node();
    const zoomLayer = svg.append<SVGGElement>('g').attr('class', 'zoom-layer').style('will-change', 'transform');
    return { svg, zoomLayer };
  }

  // ── Draw all layers ──

  function drawAll(
    zoomLayer: d3.Selection<SVGGElement, unknown, null, undefined>,
    data: D3FigureData,
    svgW: number
  ): void {
    const gComp  = zoomLayer.append<SVGGElement>('g').attr('class', 'g-comp');
    const gGuide = zoomLayer.append<SVGGElement>('g').attr('class', 'g-guide');
    const gInner = zoomLayer.append<SVGGElement>('g').attr('class', 'g-inner').style('display', 'none');
    const gRail  = zoomLayer.append<SVGGElement>('g').attr('class', 'g-rail');
    const gRHub  = zoomLayer.append<SVGGElement>('g').attr('class', 'g-rhub');
    const gRoute = zoomLayer.append<SVGGElement>('g').attr('class', 'g-route');
    const gStn   = zoomLayer.append<SVGGElement>('g').attr('class', 'g-stn');
    const gBHub  = zoomLayer.append<SVGGElement>('g').attr('class', 'g-bhub');
    // g-met is inside the zoom layer so it tracks pan/zoom correctly, but
    // sits below g-rxn so reactions always paint on top of metabolite pills.
    gMetSel      = zoomLayer.append<SVGGElement>('g').attr('class', 'g-met');
    const gRxn   = zoomLayer.append<SVGGElement>('g').attr('class', 'g-rxn');
    const gLbl   = zoomLayer.append<SVGGElement>('g').attr('class', 'g-lbl').style('display', 'none');
    const gHL    = zoomLayer.append<SVGGElement>('g').attr('class', 'g-hl');

    const view0 = data.views[0];

    // Compartment hulls
    gComp.selectAll<SVGPathElement, Compartment>('path')
      .data(data.compartments).join('path')
      .attr('d', c => hullPath(c.hull_vertices, xS, yS))
      .attr('fill', c => c.fill)
      .attr('stroke', c => c.color)
      .attr('stroke-width', 1.8)
      .attr('stroke-dasharray', '5,3')
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'none');

    gComp.selectAll<SVGTextElement, Compartment>('text')
      .data(data.compartments).join('text')
      .attr('x', c => xS(c.label_x))
      .attr('y', c => yS(c.label_y))
      .attr('text-anchor', 'middle')
      .attr('dominant-baseline', 'auto')
      .attr('fill', c => c.color)
      .attr('font-size', 12)
      .attr('font-weight', 'bold')
      .attr('font-family', "'Inter',Arial,sans-serif")
      .style('pointer-events', 'none')
      .text(c => c.label);

    gGuide.selectAll<SVGPathElement, GuideLink>('path')
      .data(data.guide_links).join('path')
      .attr('d', g => `M${xS(g.x[0]).toFixed(1)} ${yS(g.y[0]).toFixed(1)} L${xS(g.x[1]).toFixed(1)} ${yS(g.y[1]).toFixed(1)}`)
      .attr('stroke', 'rgba(139,148,158,0.55)')
      .attr('stroke-width', 0.9)
      .attr('fill', 'none')
      .attr('stroke-dasharray', '4,3')
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'none');

    const railStyleById = new Map<string, RailRouteView>();
    view0.rail_routes.forEach(rt => railStyleById.set(rt.line_id, rt));

    gRail.selectAll<SVGPathElement, RailRouteGeom>('path')
      .data(data.rail_routes).join('path')
      .attr('class', 'rail-line')
      .attr('d', d => railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k))
      .attr('stroke', d => railStyleById.get(d.line_id)?.color ?? 'rgba(110,118,129,0.18)')
      .attr('stroke-width', d => railStyleById.get(d.line_id)?.width ?? 0.8)
      .attr('fill', 'none')
      .attr('stroke-linecap', 'round')
      .attr('stroke-linejoin', 'round')
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'none');

    gMetHitSel = gMetSel!.selectAll<SVGPathElement, MetNode>('path.met-hit')
      .data(data.met_nodes).join('path')
      .attr('class', 'met-hit')
      .attr('d', d => metPath(d))
      .attr('transform', d => metTransform(d, currentTransform))
      .attr('fill', 'transparent')
      .attr('stroke', 'transparent')
      .attr('stroke-width', MET_HOVER_STROKE)
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'all')
      .style('cursor', 'default')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, metTooltipHTML(d), d.comp_color))
      .on('mousemove', moveTT)
      .on('mouseout', hideTT);

    gMetSel!.selectAll<SVGPathElement, MetNode>('path.met-shape')
      .data(data.met_nodes).join('path')
      .attr('class', 'met-shape')
      .attr('d', d => metPath(d))
      .attr('transform', d => metTransform(d, currentTransform))
      .attr('fill', d => d.display_kind === 'stop' ? '#ffffff' : d.comp_color)
      .attr('fill-opacity', d => d.display_kind === 'stop' ? 1.0 : 0.85)
      .attr('stroke', d => d.display_kind === 'stop' ? d.comp_color : '#0d1117')
      .attr('stroke-width', d => d.display_kind === 'stop' ? 1.3 : 0.7)
      .style('cursor', 'default')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, metTooltipHTML(d), d.comp_color))
      .on('mousemove', moveTT)
      .on('mouseout',  hideTT);

    // ── Routes (railway lines between hub stations) ──
    // Each station emits one path per metabolite-class line. Width and color
    // come from the per-view StationViewEntry; geometry is pinned across views.
    interface RouteDatum {
      pair_id: string;
      klass: string;
      path: [number, number][];
      lineIndex: number;
    }
    const routeData: RouteDatum[] = [];
    data.stations.forEach(st => {
      st.lines.forEach((ln, li) => {
        routeData.push({
          pair_id: st.pair_id, klass: ln.klass,
          path: ln.path, lineIndex: li,
        });
      });
    });
    const stationByPairView = new Map<string, StationViewEntry>();
    view0.stations.forEach(s => stationByPairView.set(s.pair_id, s));
    function lineStyle(d: RouteDatum, view: ViewData): { color: string; width: number } {
      const sv = view.stations.find(s => s.pair_id === d.pair_id);
      if (!sv || !sv.lines[d.lineIndex]) {
        return { color: '#6e7681', width: 2 };
      }
      const lv = sv.lines[d.lineIndex];
      return { color: lv.color, width: lv.width };
    }

    gRoute.selectAll<SVGPathElement, RouteDatum>('path')
      .data(routeData).join('path')
      .attr('class', 'route-line')
      .attr('d', d => polylinePath(d.path, xS, yS))
      .attr('fill', 'none')
      .attr('stroke-linejoin', 'round')
      .attr('stroke-linecap', 'round')
      .attr('vector-effect', 'non-scaling-stroke')
      .attr('stroke', d => lineStyle(d, view0).color)
      .attr('stroke-width', d => lineStyle(d, view0).width)
      .style('pointer-events', 'none');

    // ── Stations (vertical pills at compartment hull edges) ──
    interface PillDatum { st: StationGeom; side: 'a' | 'b'; }
    const pillData: PillDatum[] = [];
    data.stations.forEach(st => {
      pillData.push({ st, side: 'a' });
      pillData.push({ st, side: 'b' });
    });
    function pillRotation(p: PillDatum): number {
      // Pill long axis = aligned with the local hull tangent (i.e. parallel
      // to the compartment edge), so the rail line meets the pill broadside —
      // matches Birmingham-tube station pill orientation.
      const t = p.side === 'a' ? p.st.tangent_a : p.st.tangent_b;
      // Screen y is inverted (yS scale flips), so negate ty for screen-space angle.
      // +90° puts the pill's long (vertical-in-local) axis along the tangent.
      const ang = Math.atan2(-t[1], t[0]) * 180 / Math.PI + 90;
      return ang;
    }
    function pillTransform(p: PillDatum): string {
      const a = p.side === 'a' ? p.st.anchor_a : p.st.anchor_b;
      return `translate(${xS(a[0]).toFixed(1)},${yS(a[1]).toFixed(1)}) rotate(${pillRotation(p).toFixed(1)})`;
    }
    // Pill: drawn centred at origin; long axis = pill_height (vertical),
    // short axis = pill_width (horizontal). The user-facing requirement is
    // "larger in the vertical direction than in the horizontal direction".
    // Use a fixed pixel scale so pills don't shrink to invisibility on small
    // grids; convert model-space pill_height to pixels via the scale.
    function pillPixelHeight(p: PillDatum): number {
      const px0 = yS(0), px1 = yS(p.st.pill_height);
      return Math.max(7, Math.abs(px1 - px0));
    }
    function pillPixelWidth(p: PillDatum): number {
      return Math.max(3, pillPixelHeight(p) * 0.35);
    }

    gStn.selectAll<SVGRectElement, PillDatum>('rect.stn-pill')
      .data(pillData).join('rect')
      .attr('class', 'stn-pill')
      .attr('x', p => -pillPixelWidth(p) / 2)
      .attr('y', p => -pillPixelHeight(p) / 2)
      .attr('width',  p => pillPixelWidth(p))
      .attr('height', p => pillPixelHeight(p))
      .attr('rx', p => pillPixelWidth(p) / 2)
      .attr('ry', p => pillPixelWidth(p) / 2)
      .attr('fill', '#ffffff')
      .attr('stroke', p => p.side === 'a' ? p.st.comp_a_color : p.st.comp_b_color)
      .attr('stroke-width', 2.2)
      .attr('vector-effect', 'non-scaling-stroke')
      .attr('transform', pillTransform)
      .style('cursor', 'pointer')
      .on('mouseover', (ev: MouseEvent, p) =>
        showTT(ev, stationTooltipHTML(p.st),
          p.side === 'a' ? p.st.comp_a_color : p.st.comp_b_color))
      .on('mousemove', moveTT)
      .on('mouseout', hideTT);
    // ── Branch hubs (small pills at divergence points) ──
    gBHub.selectAll<SVGRectElement, BranchHub>('rect.bhub')
      .data(data.branch_hubs).join('rect')
      .attr('class', 'bhub')
      .attr('x', -4).attr('y', -10)
      .attr('width', 8).attr('height', 20)
      .attr('rx', 4).attr('ry', 4)
      .attr('fill', '#ffffff')
      .attr('stroke', '#8b949e')
      .attr('stroke-width', 1.6)
      .attr('vector-effect', 'non-scaling-stroke')
      .attr('transform', b => `translate(${xS(b.x).toFixed(1)},${yS(b.y).toFixed(1)})`)
      .style('cursor', 'pointer')
      .on('mouseover', (ev: MouseEvent, b) => showTT(ev, branchTooltipHTML(b), '#8b949e'))
      .on('mousemove', moveTT)
      .on('mouseout', hideTT);

    // ── Inner edges (transport/exchange rxn ↔ metabolite, hidden by default) ──
    gInner.selectAll<SVGPathElement, InnerEdge>('path')
      .data(data.inner_edges).join('path')
      .attr('d', e => `M${xS(e.x[0]).toFixed(1)} ${yS(e.y[0]).toFixed(1)} L${xS(e.x[1]).toFixed(1)} ${yS(e.y[1]).toFixed(1)}`)
      .attr('stroke', 'rgba(139,148,158,0.50)')
      .attr('stroke-width', 0.8)
      .attr('stroke-dasharray', '3,3')
      .attr('fill', 'none')
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'none');

    // Reactions
    gRxnHitSel = gRxn.selectAll<SVGPathElement, RxnNode>('path.rxn-hit')
      .data(view0.rxn_nodes).join('path')
      .attr('class', 'rxn-hit')
      .attr('d', d => symPath(d.symbol, d.r))
      .attr('transform', d => rxnTransform(d))
      .attr('fill', 'transparent')
      .attr('stroke', 'transparent')
      .attr('stroke-width', RXN_HOVER_STROKE)
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'all')
      .style('cursor', 'pointer')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, rxnTooltipHTML(d.hover), d.border_color))
      .on('mousemove', moveTT)
      .on('mouseout', hideTT);

    gRxnSel = gRxn.selectAll<SVGPathElement, RxnNode>('path.rxn-shape')
      .data(view0.rxn_nodes).join('path')
      .attr('class', 'rxn-shape')
      .attr('d', d => symPath(d.symbol, d.r))
      .attr('transform', d => rxnTransform(d))
      .attr('fill', d => d.fill_color)
      .attr('stroke', d => d.border_color)
      .attr('stroke-width', d => d.border_width)
      .attr('opacity', d => d.opacity)
      .style('cursor', 'pointer')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, rxnTooltipHTML(d.hover), d.border_color))
      .on('mousemove', moveTT)
      .on('mouseout',  hideTT);

    // Labels
    gLbl.selectAll<SVGTextElement, RxnNode>('text')
      .data(view0.rxn_nodes).join('text')
      .attr('x', d => xS(d.x))
      .attr('y', d => yS(d.y) - d.r - 2)
      .attr('text-anchor', 'middle')
      .attr('fill', '#8b949e')
      .attr('font-size', 7)
      .attr('font-family', "'Inter',Arial,sans-serif")
      .style('pointer-events', 'none')
      .text(d => d.name);

    // Highlight circle (search)
    gHL.append('circle')
      .attr('class', 'flov-hl')
      .attr('r', 14)
      .attr('fill', 'none')
      .attr('stroke', '#f0a030')
      .attr('stroke-width', 3)
      .attr('opacity', 0)
      .style('pointer-events', 'none');
  }

  function rxnTransform(d: RxnNode, transform: d3.ZoomTransform = currentTransform): string {
    const tx = xS(d.x).toFixed(1), ty = yS(d.y).toFixed(1);
    const inv = (1 / transform.k).toFixed(6);
    return d.symbol === 'triangle-down'
      ? `matrix(-${inv},0,0,-${inv},${tx},${ty})`
      : `matrix(${inv},0,0,${inv},${tx},${ty})`;
  }

  // ── Toolbar ──

  function buildToolbar(data: D3FigureData, parent: HTMLElement): void {
    const bar = document.createElement('div');
    bar.className = 'flov-toolbar';
    parent.appendChild(bar);

    // Flux view group
    if (data.view_labels.length) {
      const grp = btnGroup(bar, 'Flux view');
      data.view_labels.forEach((lbl, i) => {
        const btn = mkBtn(lbl, i === 0);
        btn.onclick = () => switchView(i);
        grp.appendChild(btn);
        viewBtns.push(btn);
      });
    }

    // Metabolites
    mkToggle(bar, 'Metabolites', 'Show', 'Hide', true, on => {
      gMetSel!.style('display', on ? '' : 'none');
    });

    // Inner connections (rxn → metabolite edges of transport/exchange reactions)
    mkToggle(bar, 'Inner connections', 'Show', 'Hide', false, on => {
      const el = wrapper.querySelector<HTMLElement>('.g-inner');
      if (el) el.style.display = on ? '' : 'none';
    });

    // Dark/light mode
    const themeGroup = btnGroup(bar, 'Theme');
    const themeBtn = mkBtn('Dark mode', true);
    const applyTheme = (): void => {
      wrapper.classList.toggle('light-mode', lightMode);
      if (svgEl) svgEl.style.background = lightMode ? '#f6f8fa' : '#0d1117';
      themeBtn.textContent = lightMode ? 'Dark mode' : 'Light mode';
      themeBtn.classList.toggle('active', lightMode);
    };
    themeBtn.onclick = () => {
      lightMode = !lightMode;
      applyTheme();
    };
    applyTheme();
    themeGroup.appendChild(themeBtn);

    // Reset zoom
    const rg = btnGroup(bar, '');
    const rb = mkBtn('⌖ Reset zoom', false);
    rb.onclick = () => {
      if (svgEl && zoomB)
        d3.select(svgEl).transition().duration(350).call(zoomB.transform, d3.zoomIdentity);
    };
    rg.appendChild(rb);
  }

  function btnGroup(bar: HTMLElement, label: string): HTMLElement {
    const grp = document.createElement('div');
    grp.className = 'flov-btn-group';
    if (label) {
      const lbl = document.createElement('span');
      lbl.className = 'flov-btn-label';
      lbl.textContent = label;
      grp.appendChild(lbl);
    }
    bar.appendChild(grp);
    return grp;
  }

  function mkBtn(text: string, active: boolean): HTMLButtonElement {
    const b = document.createElement('button');
    b.className = 'flov-btn' + (active ? ' active' : '');
    b.textContent = text;
    return b;
  }

  function mkToggle(
    bar: HTMLElement, label: string,
    onTxt: string, offTxt: string,
    initOn: boolean,
    cb: (on: boolean) => void
  ): void {
    const grp = btnGroup(bar, label);
    const bOn  = mkBtn(onTxt,  initOn);
    const bOff = mkBtn(offTxt, !initOn);
    bOn.onclick  = () => { bOn.classList.add('active');  bOff.classList.remove('active'); cb(true); };
    bOff.onclick = () => { bOff.classList.add('active'); bOn.classList.remove('active');  cb(false); };
    grp.appendChild(bOn);
    grp.appendChild(bOff);
  }

  // ── View switching ──

  function switchView(vi: number): void {
    if (!fig) return;
    currentView = vi;
    const view = fig.views[vi];
    const zl = wrapper.querySelector<SVGGElement>('.zoom-layer');
    if (!zl) return;

    const railViewById = new Map<string, RailRouteView>();
    view.rail_routes.forEach(rt => railViewById.set(rt.line_id, rt));
    d3.select(zl).select('.g-rail').selectAll<SVGPathElement, RailRouteGeom>('path.rail-line')
      .each(function(d: RailRouteGeom) {
        const rv = railViewById.get(d.line_id);
        d3.select(this)
          .attr('stroke', rv?.color ?? 'rgba(110,118,129,0.18)')
          .attr('stroke-width', rv?.width ?? 0.8);
      });

    d3.select(zl).select('.g-rxn').selectAll<SVGPathElement, RxnNode>('path.rxn-shape')
      .data(view.rxn_nodes)
      .attr('d', d => symPath(d.symbol, d.r))
      .attr('transform', d => rxnTransform(d))
      .attr('fill', d => d.fill_color)
      .attr('opacity', d => d.opacity);

    // Update route line widths & colors for the new view
    d3.select(zl).select('.g-route').selectAll<SVGPathElement, {pair_id:string; lineIndex:number}>('path.route-line')
      .each(function(d: {pair_id:string; lineIndex:number}) {
        const sv = view.stations.find(s => s.pair_id === d.pair_id);
        const lv = sv ? sv.lines[d.lineIndex] : undefined;
        if (lv) {
          d3.select(this).attr('stroke', lv.color).attr('stroke-width', lv.width);
        } else {
          d3.select(this).attr('stroke', '#6e7681').attr('stroke-width', 1);
        }
      });

    gRxnHitSel!
      .data(view.rxn_nodes)
      .attr('d', d => symPath(d.symbol, d.r))
      .attr('transform', d => rxnTransform(d));

    d3.select(zl).select('.g-lbl').selectAll<SVGTextElement, RxnNode>('text')
      .data(view.rxn_nodes)
      .attr('x', d => xS(d.x))
      .attr('y', d => yS(d.y) - d.r - 2);

    viewBtns.forEach((b, i) => b.classList.toggle('active', i === vi));
  }

  // ── Tooltip ──

  function buildTooltip(parent: HTMLElement): void {
    const el = document.createElement('div');
    el.className = 'flov-tt';
    parent.appendChild(el);
    tooltipEl = el;
  }

  function showTT(ev: MouseEvent, html: string, borderColor: string): void {
    if (!tooltipEl) return;
    tooltipEl.innerHTML = html;
    tooltipEl.style.borderColor = borderColor;
    tooltipEl.style.display = 'block';
    moveTT(ev);
  }

  function moveTT(ev: MouseEvent): void {
    if (!tooltipEl) return;
    // Always read a fresh rect — caching across resizes / scrolls / dev-tool
    // toggles puts the tooltip far from the cursor. getBoundingClientRect is
    // cheap (single DOM read) and runs only on mousemove inside the SVG.
    const rect = wrapper.getBoundingClientRect();
    let tx = ev.clientX - rect.left + 14;
    let ty = ev.clientY - rect.top - 10;
    const tw = tooltipEl.offsetWidth || 220;
    const th = tooltipEl.offsetHeight || 110;
    if (tx + tw > rect.width - 6)  tx = ev.clientX - rect.left - tw - 14;
    if (ty + th > rect.height - 6) ty = ev.clientY - rect.top - th - 10;
    tooltipEl.style.left = tx + 'px';
    tooltipEl.style.top  = ty + 'px';
  }

  function hideTT(): void {
    if (tooltipEl) tooltipEl.style.display = 'none';
  }

  // ── Colorbar ──

  function buildColorbar(parent: HTMLElement, absMax: number, svgH: number): void {
    const cbH = 210, cbW = 52;
    const cb = d3.select(parent).append<SVGSVGElement>('svg')
      .attr('width', cbW).attr('height', cbH)
      .style('position', 'absolute')
      .style('right', '8px')
      .style('top', Math.max(52, (svgH - cbH) / 2 + 40) + 'px');

    const defs = cb.append('defs');
    const grad = defs.append('linearGradient')
      .attr('id', 'flov-cb-grad')
      .attr('x1', '0').attr('y1', '0').attr('x2', '0').attr('y2', '1');
    // top → bottom: dark red → light red → grey → light blue → dark blue
    [
      ['0%',   '#c62828'],
      ['25%',  '#ef9a9a'],
      ['50%',  '#6e7681'],
      ['75%',  '#90caf9'],
      ['100%', '#1565c0'],
    ].forEach(([off, col]) => grad.append('stop').attr('offset', off).attr('stop-color', col));

    cb.append('rect')
      .attr('x', 2).attr('y', 12)
      .attr('width', 18).attr('height', cbH - 24)
      .attr('fill', 'url(#flov-cb-grad)')
      .attr('rx', 3);

    const sc = d3.scaleLinear().domain([absMax, -absMax]).range([12, cbH - 12]);
    cb.append<SVGGElement>('g')
      .attr('transform', 'translate(20,0)')
      .call(
        d3.axisRight(sc).ticks(5)
          .tickFormat(d3.format('+.1f') as (d: d3.NumberValue) => string)
      )
      .call(g => {
        g.select('.domain').attr('stroke', '#444c56');
        g.selectAll<SVGTextElement, unknown>('text').attr('fill', '#8b949e').attr('font-size', 9);
        g.selectAll<SVGLineElement, unknown>('.tick line').attr('stroke', '#444c56');
      });

    cb.append('text')
      .attr('x', 11).attr('y', 9)
      .attr('text-anchor', 'middle')
      .attr('fill', '#7d8590')
      .attr('font-size', 8)
      .attr('font-family', "'Inter',Arial,sans-serif")
      .text('mmol/gDW/h');
  }

  // ── Search ──

  function buildSearch(
    parent: HTMLElement, data: D3FigureData,
    svgNode: SVGSVGElement, svgW: number, svgH: number
  ): void {
    const ia = data.interactivity;
    const rxnIdsLc   = ia.rxn_ids.map(s => s.toLowerCase());
    const rxnNamesLc = ia.rxn_names.map(s => s.toLowerCase());
    const metIdsLc   = ia.met_ids.map(s => s.toLowerCase());
    const metNamesLc = ia.met_names.map(s => s.toLowerCase());

    const wrap = document.createElement('div');
    wrap.className = 'flov-search-wrap';
    const inp = document.createElement('input');
    inp.className = 'flov-search-input';
    inp.placeholder = '🔍  Search reaction / metabolite…';
    inp.type = 'text';
    const res = document.createElement('div');
    res.className = 'flov-search-results';
    wrap.appendChild(inp);
    wrap.appendChild(res);
    parent.appendChild(wrap);

    const hlCircle = (): SVGCircleElement | null =>
      wrapper.querySelector<SVGCircleElement>('.flov-hl');

    function clearHL(): void {
      const c = hlCircle(); if (c) c.setAttribute('opacity', '0');
    }

    function setHL(px: number, py: number): void {
      const c = hlCircle();
      if (!c) return;
      c.setAttribute('cx', String(px));
      c.setAttribute('cy', String(py));
      c.setAttribute('opacity', '1');
    }

    function zoomToPoint(px: number, py: number): void {
      const k = 5;
      const t = d3.zoomIdentity.translate(svgW / 2 - px * k, svgH / 2 - py * k).scale(k);
      d3.select(svgNode).transition().duration(400).call(zoomB!.transform, t);
    }

    function zoomToFit(pxs: number[], pys: number[]): void {
      const x0 = Math.min(...pxs), x1 = Math.max(...pxs);
      const y0 = Math.min(...pys), y1 = Math.max(...pys);
      const pad = 60;
      const k = Math.min(8, svgW / (x1 - x0 + pad * 2), svgH / (y1 - y0 + pad * 2));
      const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
      const t = d3.zoomIdentity.translate(svgW / 2 - cx * k, svgH / 2 - cy * k).scale(k);
      d3.select(svgNode).transition().duration(400).call(zoomB!.transform, t);
    }

    function rxnPX(i: number): number { return xS(ia.rxn_x[i]); }
    function rxnPY(i: number): number { return yS(ia.rxn_y[i]); }
    function metPX(i: number): number { return xS(ia.met_x[i]); }
    function metPY(i: number): number { return yS(ia.met_y[i]); }

    function selectResult(type: 'rxn' | 'met', idx: number): void {
      const px = type === 'rxn' ? rxnPX(idx) : metPX(idx);
      const py = type === 'rxn' ? rxnPY(idx) : metPY(idx);
      setHL(px, py);
      inp.value = type === 'rxn'
        ? (ia.rxn_names[idx] || ia.rxn_ids[idx])
        : (ia.met_names[idx] || ia.met_ids[idx]);
      res.style.display = 'none';
      zoomToPoint(px, py);
    }

    function runSearch(q: string): void {
      res.innerHTML = '';
      if (!q.trim()) { clearHL(); res.style.display = 'none'; return; }
      const lq = q.toLowerCase();
      const hpxs: number[] = [], hpys: number[] = [];
      const hits: { type: 'rxn' | 'met'; idx: number }[] = [];

      for (let i = 0; i < ia.rxn_ids.length; i++) {
        if (rxnIdsLc[i].includes(lq) || rxnNamesLc[i].includes(lq)) {
          hpxs.push(rxnPX(i)); hpys.push(rxnPY(i));
          if (hits.length < 20) hits.push({ type: 'rxn', idx: i });
        }
      }
      for (let i = 0; i < ia.met_ids.length; i++) {
        if (metIdsLc[i].includes(lq) || metNamesLc[i].includes(lq)) {
          hpxs.push(metPX(i)); hpys.push(metPY(i));
          if (hits.length < 20) hits.push({ type: 'met', idx: i });
        }
      }

      if (!hpxs.length) {
        const d = document.createElement('div');
        d.className = 'flov-search-result';
        d.style.color = '#7d8590';
        d.textContent = 'No matches';
        res.appendChild(d);
        res.style.display = 'block';
        clearHL();
        return;
      }

      if (hpxs.length === 1) zoomToPoint(hpxs[0], hpys[0]);
      else zoomToFit(hpxs, hpys);

      if (hits.length === 1) {
        const h = hits[0];
        setHL(h.type === 'rxn' ? rxnPX(h.idx) : metPX(h.idx),
               h.type === 'rxn' ? rxnPY(h.idx) : metPY(h.idx));
      } else {
        clearHL();
      }

      for (const h of hits) {
        const id   = h.type === 'rxn' ? ia.rxn_ids[h.idx]   : ia.met_ids[h.idx];
        const name = h.type === 'rxn' ? ia.rxn_names[h.idx]  : ia.met_names[h.idx];
        const el = document.createElement('div');
        el.className = 'flov-search-result';
        el.title = id + (name && name !== id ? ' — ' + name : '');
        el.innerHTML =
          `<b>${id}</b>` +
          (name && name !== id ? ` <span style="color:#7d8590">— ${name}</span>` : '') +
          ` <span style="color:#444c56;font-size:10px">[${h.type}]</span>`;
        el.onclick = () => selectResult(h.type, h.idx);
        res.appendChild(el);
      }
      res.style.display = 'block';
    }

    inp.addEventListener('input', function (this: HTMLInputElement) { runSearch(this.value); });
    inp.addEventListener('keydown', (e: KeyboardEvent) => {
      if (e.key === 'Escape') { inp.value = ''; clearHL(); res.style.display = 'none'; }
    });
    // Self-removing click-outside handler
    document.addEventListener('click', function outside(e: MouseEvent) {
      if (!wrap.isConnected) { document.removeEventListener('click', outside); return; }
      if (!wrap.contains(e.target as Node)) res.style.display = 'none';
    });
  }

  // ── Resize handler ──

  new ResizeObserver((entries: ResizeObserverEntry[]) => {
    for (const entry of entries) {
      const w = entry.contentRect.width;
      if (w <= 0) continue;
      cachedWrapperRect = null;
      if (!fig) continue;
      clearTimeout(resizeTimer!);
      resizeTimer = setTimeout(() => buildUI(fig!), 120);
    }
  }).observe(I);

  // ── Model listener ──

  function onData(): void {
    const raw: string = S.get('_figure_json');
    if (!raw || raw === '{}') return;
    let data: D3FigureData;
    try { data = JSON.parse(raw); } catch { return; }
    buildUI(data);
  }

  S.on('change:_figure_json', onData);
  onData();
}

export default { render };
