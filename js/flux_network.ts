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
interface MetPosition { x: number; y: number; }
interface EdgeGroup {
  class: 'reg' | 'ss';
  x: (number | null)[]; y: (number | null)[];
  color: string; width: number;
}
interface ViewData {
  label: string; rxn_nodes: RxnNode[];
  met_positions: MetPosition[]; edge_groups: EdgeGroup[];
}
interface MetDatum {
  node: MetNode; pos: MetPosition;
}
interface MetNode {
  id: string; name: string;
  comp: string; comp_label: string; comp_color: string;
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
}
interface Interactivity {
  view_labels: string[];
  rxn_ids: string[]; rxn_names: string[]; rxn_x: number[]; rxn_y: number[];
  met_ids: string[]; met_names: string[];
  met_flux_x: number[][]; met_flux_y: number[][];
}
interface D3FigureData {
  meta: Meta; compartments: Compartment[];
  view_labels: string[]; views: ViewData[];
  met_nodes: MetNode[]; interactivity: Interactivity;
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

// Fixed 4-pixel-radius diamond for metabolites
function metPath(): string {
  return d3.symbol(d3.symbolDiamond, Math.PI * 16)(null) ?? '';
}

// ── Path utilities ────────────────────────────────────────────────────────────

function edgePath(
  xs: (number | null)[], ys: (number | null)[],
  xS: d3.ScaleLinear<number, number>, yS: d3.ScaleLinear<number, number>
): string {
  let d = '', open = false;
  for (let i = 0; i < xs.length; i++) {
    if (xs[i] === null) { open = false; continue; }
    const px = xS(xs[i]!).toFixed(1), py = yS(ys[i]!).toFixed(1);
    d += open ? ` L${px} ${py}` : `M${px} ${py}`;
    open = true;
  }
  return d;
}

function hullPath(
  verts: [number, number][],
  xS: d3.ScaleLinear<number, number>, yS: d3.ScaleLinear<number, number>
): string {
  const gen = d3.line<[number, number]>()
    .x(d => xS(d[0])).y(d => yS(d[1]))
    .curve(d3.curveCatmullRomClosed.alpha(0.5));
  return gen(verts) ?? '';
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
  let gMetHitSel: d3.Selection<SVGPathElement, MetDatum, SVGGElement, unknown> | null = null;
  let gRxnSel: d3.Selection<SVGPathElement, RxnNode, SVGGElement, unknown> | null = null;
  let gRxnHitSel: d3.Selection<SVGPathElement, RxnNode, SVGGElement, unknown> | null = null;
  let currentTransform: d3.ZoomTransform = d3.zoomIdentity;
  let rafId: number | null = null;
  let cachedWrapperRect: DOMRect | null = null;
  let lightMode = true;

  const MARGIN = { t: 10, r: 90, b: 10, l: 20 };
  const MET_HOVER_STROKE = 12;
  const RXN_HOVER_STROKE = 14;

  function metTransform(pos: MetPosition, transform: d3.ZoomTransform): string {
    const sx = transform.applyX(xS(pos.x)).toFixed(1);
    const sy = transform.applyY(yS(pos.y)).toFixed(1);
    return `translate(${sx},${sy})`;
  }

  function updateMetZoom(transform: d3.ZoomTransform): void {
    if (!gMetSel) return;
    gMetSel.selectAll<SVGPathElement, MetDatum>('path.met-shape')
      .attr('transform', d => metTransform(d.pos, transform));
    gMetHitSel?.attr('transform', d => metTransform(d.pos, transform));
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
    gMetSel = svg.insert<SVGGElement>('g', '.zoom-layer').attr('class', 'g-met').style('display', 'none');
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
    const gSSEdge = zoomLayer.append<SVGGElement>('g').attr('class', 'g-ss').style('display', 'none');
    const gEdge  = zoomLayer.append<SVGGElement>('g').attr('class', 'g-edge');
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

    // Edges (view 0)
    const regEgs = view0.edge_groups.filter(e => e.class === 'reg');
    const ssEgs  = view0.edge_groups.filter(e => e.class === 'ss');

    gEdge.selectAll<SVGPathElement, EdgeGroup>('path')
      .data(regEgs).join('path')
      .attr('d', eg => edgePath(eg.x, eg.y, xS, yS))
      .attr('stroke', eg => eg.color)
      .attr('stroke-width', eg => eg.width)
      .attr('fill', 'none')
      .attr('stroke-linecap', 'round')
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'none');

    gSSEdge.selectAll<SVGPathElement, EdgeGroup>('path')
      .data(ssEgs).join('path')
      .attr('d', eg => edgePath(eg.x, eg.y, xS, yS))
      .attr('stroke', eg => eg.color)
      .attr('stroke-width', eg => eg.width)
      .attr('fill', 'none')
      .attr('stroke-dasharray', '4,3')
      .attr('stroke-linecap', 'round')
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'none');

    // Metabolites (initially hidden)
    const metData: MetDatum[] = data.met_nodes.map((node, i) => ({
      node, pos: view0.met_positions[i],
    }));

    gMetHitSel = gMetSel!.selectAll<SVGPathElement, MetDatum>('path.met-hit')
      .data(metData).join('path')
      .attr('class', 'met-hit')
      .attr('d', metPath())
      .attr('transform', d => metTransform(d.pos, currentTransform))
      .attr('fill', 'transparent')
      .attr('stroke', 'transparent')
      .attr('stroke-width', MET_HOVER_STROKE)
      .attr('vector-effect', 'non-scaling-stroke')
      .style('pointer-events', 'all')
      .style('cursor', 'default')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, metTooltipHTML(d.node), d.node.comp_color))
      .on('mousemove', moveTT)
      .on('mouseout', hideTT);

    gMetSel!.selectAll<SVGPathElement, MetDatum>('path.met-shape')
      .data(metData).join('path')
      .attr('class', 'met-shape')
      .attr('d', metPath())
      .attr('transform', d => metTransform(d.pos, currentTransform))
      .attr('fill', d => d.node.comp_color)
      .attr('fill-opacity', 0.85)
      .attr('stroke', '#0d1117')
      .attr('stroke-width', 0.7)
      .style('cursor', 'default')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, metTooltipHTML(d.node), d.node.comp_color))
      .on('mousemove', moveTT)
      .on('mouseout',  hideTT);

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
    mkToggle(bar, 'Metabolites', 'Show', 'Hide', false, on => {
      gMetSel!.style('display', on ? '' : 'none');
    });

    // Exchange / Transp.
    mkToggle(bar, 'Exchange / Transp.', 'Show', 'Hide', false, on => {
      wrapper.querySelector<HTMLElement>('.g-ss')!.style.display = on ? '' : 'none';
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

    const regEgs = view.edge_groups.filter(e => e.class === 'reg');
    const ssEgs  = view.edge_groups.filter(e => e.class === 'ss');

    d3.select(zl).select('.g-edge').selectAll<SVGPathElement, EdgeGroup>('path')
      .data(regEgs)
      .attr('d', eg => edgePath(eg.x, eg.y, xS, yS))
      .attr('stroke', eg => eg.color)
      .attr('stroke-width', eg => eg.width);

    d3.select(zl).select('.g-ss').selectAll<SVGPathElement, EdgeGroup>('path')
      .data(ssEgs)
      .attr('d', eg => edgePath(eg.x, eg.y, xS, yS))
      .attr('stroke', eg => eg.color)
      .attr('stroke-width', eg => eg.width);

    d3.select(zl).select('.g-rxn').selectAll<SVGPathElement, RxnNode>('path.rxn-shape')
      .data(view.rxn_nodes)
      .attr('d', d => symPath(d.symbol, d.r))
      .attr('transform', d => rxnTransform(d))
      .attr('fill', d => d.fill_color)
      .attr('opacity', d => d.opacity);

    const metData: MetDatum[] = fig.met_nodes.map((node, i) => ({
      node, pos: view.met_positions[i],
    }));
    gMetSel!.selectAll<SVGPathElement, MetDatum>('path.met-shape')
      .data(metData)
      .attr('transform', d => metTransform(d.pos, currentTransform));
    gMetHitSel!
      .data(metData)
      .attr('transform', d => metTransform(d.pos, currentTransform));

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
    if (!cachedWrapperRect) cachedWrapperRect = wrapper.getBoundingClientRect();
    const rect = cachedWrapperRect;
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
    function metPX(i: number): number { return xS(ia.met_flux_x[currentView][i]); }
    function metPY(i: number): number { return yS(ia.met_flux_y[currentView][i]); }

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
