import * as d3 from 'd3';
import type {
  Compartment,
  D3FigureData,
  FocusState,
  GuideLink,
  InnerEdge,
  MetNode,
  RailRouteGeom,
  RailRouteView,
  RouteDatum,
  RxnNode,
  RxnSummary,
  StationGeom,
  StationLineGeom,
  ViewData,
} from './types';
import {
  injectCSS,
  hullPath,
  KLASS_COLORS,
  metPath,
  metTooltipHTML,
  railOffsetPath,
  rxnTooltipHTML,
  STATION_BRANCH_COLORS,
  symPath,
} from './rendering';

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
  let stationPanelEl: HTMLDivElement | null = null;
  let activeStationPanel: StationGeom | null = null;
  let currentStationSearch = '';
  let resizeTimer: ReturnType<typeof setTimeout> | null = null;
  let gMetSel: d3.Selection<SVGGElement, unknown, null, undefined> | null = null;
  let gMetHitSel: d3.Selection<SVGPathElement, MetNode, SVGGElement, unknown> | null = null;
  let gRxnSel: d3.Selection<SVGPathElement, RxnNode, SVGGElement, unknown> | null = null;
  let gRxnHitSel: d3.Selection<SVGPathElement, RxnNode, SVGGElement, unknown> | null = null;
  let currentTransform: d3.ZoomTransform = d3.zoomIdentity;
  let rafId: number | null = null;
  let lightMode = true;
  let focusState: FocusState | null = null;
  const arrowMarkerEndId = `flov-arrow-end-${Math.random().toString(36).slice(2)}`;
  const arrowMarkerStartId = `flov-arrow-start-${Math.random().toString(36).slice(2)}`;

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

  function updateRouteZoom(transform: d3.ZoomTransform): void {
    const zl = wrapper.querySelector<SVGGElement>('.zoom-layer');
    if (!zl || !fig) return;
    d3.select(zl).select('.g-route')
      .selectAll<SVGPathElement, RouteDatum>('path.route-line')
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
    focusState = null;
    activeStationPanel = null;
    if (rafId !== null) { cancelAnimationFrame(rafId); rafId = null; }
    wrapper.innerHTML = '';
    tooltipEl = null;
    stationPanelEl = null;
    viewBtns = [];

    const svgW = Math.max(wrapper.clientWidth || 900, 600);
    const svgH = data.meta.height;

    xS = d3.scaleLinear().domain(data.meta.x_range).range([MARGIN.l, svgW - MARGIN.r]);
    yS = d3.scaleLinear().domain(data.meta.y_range).range([svgH - MARGIN.b, MARGIN.t]);

    wrapper.classList.toggle('light-mode', lightMode);
    buildToolbar(data, wrapper);
    const { svg, zoomLayer } = buildSVG(svgW, svgH);
    buildLegend(wrapper, data, svgH);
    buildSearch(wrapper, data);
    buildStationPanel(wrapper);
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
          updateRouteZoom(currentTransform);
        });
      });
    svg.call(zoomB);
    svg.on('dblclick.zoom', () => {
      svg.transition().duration(380).call(zoomB!.transform, d3.zoomIdentity);
    });
    svg.on('click.focus', () => clearFocus());
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
    const defs = svg.append('defs');
    defs.append('marker')
      .attr('id', arrowMarkerEndId)
      .attr('viewBox', '0 0 8 8')
      .attr('refX', 7.2)
      .attr('refY', 4)
      .attr('markerWidth', 3.5)
      .attr('markerHeight', 3.5)
      .attr('orient', 'auto')
      .append('path')
      .attr('d', 'M0,0 L8,4 L0,8 Z')
      .attr('fill', 'context-stroke');
    defs.append('marker')
      .attr('id', arrowMarkerStartId)
      .attr('viewBox', '0 0 8 8')
      .attr('refX', 7.2)
      .attr('refY', 4)
      .attr('markerWidth', 3.5)
      .attr('markerHeight', 3.5)
      .attr('orient', 'auto-start-reverse')
      .append('path')
      .attr('d', 'M0,0 L8,4 L0,8 Z')
      .attr('fill', 'context-stroke');
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
    const gRoute = zoomLayer.append<SVGGElement>('g').attr('class', 'g-route');
    const gStn   = zoomLayer.append<SVGGElement>('g').attr('class', 'g-stn');
    // g-met is inside the zoom layer so it tracks pan/zoom correctly, but
    // sits below g-rxn so reactions always paint on top of metabolite pills.
    gMetSel      = zoomLayer.append<SVGGElement>('g').attr('class', 'g-met');
    const gRxn   = zoomLayer.append<SVGGElement>('g').attr('class', 'g-rxn');
    const gLbl   = zoomLayer.append<SVGGElement>('g').attr('class', 'g-lbl').style('display', 'none');

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
      .style('cursor', 'pointer')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, metTooltipHTML(d), d.comp_color))
      .on('mousemove', moveTT)
      .on('mouseout', hideTT)
      .on('click', (ev: MouseEvent, d) => setMetFocus(d, ev));

    gMetSel!.selectAll<SVGPathElement, MetNode>('path.met-shape')
      .data(data.met_nodes).join('path')
      .attr('class', 'met-shape')
      .attr('d', d => metPath(d))
      .attr('transform', d => metTransform(d, currentTransform))
      .attr('fill', d => d.display_kind === 'stop' ? '#ffffff' : d.comp_color)
      .attr('fill-opacity', d => d.display_kind === 'stop' ? 1.0 : 0.85)
      .attr('stroke', d => d.display_kind === 'stop' ? d.comp_color : '#0d1117')
      .attr('stroke-width', d => d.display_kind === 'stop' ? 1.3 : 0.7)
      .style('cursor', 'pointer')
      .on('mouseover', (ev: MouseEvent, d) => showTT(ev, metTooltipHTML(d), d.comp_color))
      .on('mousemove', moveTT)
      .on('mouseout',  hideTT)
      .on('click', (ev: MouseEvent, d) => setMetFocus(d, ev));

    // ── Routes (railway lines between hub stations) ──
    // Each station emits one path per metabolite-class line. Width and color
    // come from the per-view StationViewEntry; geometry is pinned across views.
    const routeData: RouteDatum[] = [];
    data.stations.forEach(st => {
      st.lines.forEach((ln, li) => {
        routeData.push({
          pair_id: st.pair_id, klass: ln.klass,
          rxn_ids: ln.rxn_ids,
          seg_slots: ln.seg_slots,
          path: ln.path, lineIndex: li,
        });
      });
    });
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
      .attr('d', d => railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k))
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
      .on('click', (ev: MouseEvent, p) => {
        setStationFocus(p.st, ev);
        showStationPanel(p.st);
      });
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
      .on('mouseout', hideTT)
      .on('click', (ev: MouseEvent, d) => setRxnFocus(d, ev));

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
      .on('mouseout',  hideTT)
      .on('click', (ev: MouseEvent, d) => setRxnFocus(d, ev));

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

    applyFocus();
  }

  function rxnTransform(d: RxnNode, transform: d3.ZoomTransform = currentTransform): string {
    const tx = xS(d.x).toFixed(1), ty = yS(d.y).toFixed(1);
    const inv = (1 / transform.k).toFixed(6);
    return d.symbol === 'triangle-down'
      ? `matrix(-${inv},0,0,-${inv},${tx},${ty})`
      : `matrix(${inv},0,0,${inv},${tx},${ty})`;
  }

  // ── Click focus ──

  function connectedMetIdsForRxn(rxnId: string): Set<string> {
    const ids = new Set<string>();
    if (!fig) return ids;
    fig.rail_routes.forEach(rt => { if (rt.rxn_id === rxnId) ids.add(rt.met_id); });
    fig.inner_edges.forEach(e => { if (e.rxn_id === rxnId) ids.add(e.met_id); });
    fig.met_nodes.forEach(m => {
      if (m.consumers.includes(rxnId) || m.producers.includes(rxnId)) ids.add(m.id);
    });
    return ids;
  }

  function connectedRxnIdsForMet(metId: string): Set<string> {
    const ids = new Set<string>();
    if (!fig) return ids;
    fig.rail_routes.forEach(rt => { if (rt.met_id === metId) ids.add(rt.rxn_id); });
    fig.inner_edges.forEach(e => { if (e.met_id === metId) ids.add(e.rxn_id); });
    const met = fig.met_nodes.find(m => m.id === metId);
    met?.consumers.forEach(r => ids.add(r));
    met?.producers.forEach(r => ids.add(r));
    return ids;
  }

  function reactionFocusSet(rxnId: string): { rxnIds: Set<string>; metIds: Set<string> } {
    return {
      rxnIds: new Set<string>([rxnId]),
      metIds: connectedMetIdsForRxn(rxnId),
    };
  }

  function setRxnFocus(d: RxnNode, ev: MouseEvent): void {
    ev.stopPropagation();
    if (focusState?.type === 'rxn' && focusState.id === d.id) {
      clearFocus();
      return;
    }
    const neighborhood = reactionFocusSet(d.id);
    focusState = {
      type: 'rxn',
      id: d.id,
      rxnIds: neighborhood.rxnIds,
      metIds: neighborhood.metIds,
    };
    applyFocus();
  }

  function setMetFocus(d: MetNode, ev: MouseEvent): void {
    ev.stopPropagation();
    if (focusState?.type === 'met' && focusState.id === d.id) {
      clearFocus();
      return;
    }
    focusState = {
      type: 'met',
      id: d.id,
      rxnIds: connectedRxnIdsForMet(d.id),
      metIds: new Set([d.id]),
    };
    applyFocus();
  }

  function focusMetaboliteById(metId: string): void {
    focusState = {
      type: 'met',
      id: metId,
      rxnIds: connectedRxnIdsForMet(metId),
      metIds: new Set([metId]),
    };
    applyFocus();
  }

  function setStationFocus(st: StationGeom, ev: MouseEvent): void {
    ev.stopPropagation();
    hideTT();
    const rxnIds = new Set<string>();
    const metIds = new Set<string>();
    st.lines.forEach(ln => {
      ln.rxn_ids.forEach(rid => {
        rxnIds.add(rid);
        connectedMetIdsForRxn(rid).forEach(mid => metIds.add(mid));
      });
    });
    focusState = {
      type: 'station',
      id: st.pair_id,
      rxnIds,
      metIds,
    };
    applyFocus();
  }

  function setStationReactionFocus(rxnId: string, ev: MouseEvent): void {
    ev.stopPropagation();
    focusStationReactionById(rxnId);
  }

  function focusStationReactionById(rxnId: string): void {
    const neighborhood = reactionFocusSet(rxnId);
    focusState = {
      type: 'rxn',
      id: rxnId,
      rxnIds: neighborhood.rxnIds,
      metIds: neighborhood.metIds,
    };
    applyFocus();
  }

  function clearFocus(): void {
    if (!focusState) return;
    focusState = null;
    hideStationPanel();
    applyFocus();
  }

  function focusHasStation(st: StationGeom): boolean {
    if (!focusState) return false;
    return st.lines.some(ln => ln.rxn_ids.some(rid => focusState!.rxnIds.has(rid)));
  }

  function flowDirection(
    rxnId: string,
    metId: string,
    stoich?: number
  ): 'rxn-to-met' | 'met-to-rxn' | null {
    if (typeof stoich === 'number' && Number.isFinite(stoich) && stoich !== 0) {
      return stoich > 0 ? 'rxn-to-met' : 'met-to-rxn';
    }
    const met = fig?.met_nodes.find(m => m.id === metId);
    if (met?.producers.includes(rxnId)) return 'rxn-to-met';
    if (met?.consumers.includes(rxnId)) return 'met-to-rxn';
    return null;
  }

  function flowMagnitude(rxnId: string, metId: string, stoich?: number): number {
    if (typeof stoich === 'number' && Number.isFinite(stoich) && stoich !== 0) {
      return Math.max(1, Math.abs(stoich));
    }
    if (!fig) return 1;
    let m = 0;
    for (const r of fig.rail_routes) {
      if (r.rxn_id === rxnId && r.met_id === metId) m += Math.abs(r.stoich) || 1;
    }
    if (m === 0) {
      for (const e of fig.inner_edges) {
        if (e.rxn_id === rxnId && e.met_id === metId) m += 1;
      }
    }
    return Math.max(1, m);
  }

  function applyFlowArrow(
    sel: d3.Selection<SVGPathElement, any, any, any>,
    isFocus: boolean,
    direction: 'rxn-to-met' | 'met-to-rxn' | null,
    magnitude: number = 1
  ): void {
    const activeArrow = focusState !== null && isFocus && direction !== null;
    sel.classed('flov-flow-arrow', activeArrow)
      .classed('flov-flow-forward', activeArrow && direction === 'rxn-to-met')
      .classed('flov-flow-reverse', activeArrow && direction === 'met-to-rxn')
      .style('animation-duration', activeArrow ? `${1 / magnitude}s` : '')
      .attr('marker-end', null)
      .attr('marker-start', null);
  }

  function applyFocus(): void {
    if (!fig) return;
    const zl = wrapper.querySelector<SVGGElement>('.zoom-layer');
    if (!zl) return;

    const active = focusState !== null;
    const dimStroke = lightMode ? '#8c959f' : '#484f58';
    const dimFill = lightMode ? '#d8dee4' : '#30363d';
    const view = fig.views[currentView];
    const railViewById = new Map<string, RailRouteView>();
    view.rail_routes.forEach(rt => railViewById.set(rt.line_id, rt));

    d3.select(zl).select('.g-guide').selectAll<SVGPathElement, GuideLink>('path')
      .attr('opacity', g => !active ||
        (g.node_type === 'rxn' && focusState!.rxnIds.has(g.node_id)) ||
        (g.node_type === 'met' && focusState!.metIds.has(g.node_id)) ? 1 : 0.08);

    d3.select(zl).select('.g-inner').selectAll<SVGPathElement, InnerEdge>('path')
      .each(function(e: InnerEdge) {
        const isFocus = active && focusState!.rxnIds.has(e.rxn_id) && focusState!.metIds.has(e.met_id);
        const sel = d3.select(this);
        sel
          .attr('stroke', !active || isFocus ? 'rgba(139,148,158,0.50)' : dimStroke)
          .attr('opacity', !active || isFocus ? 1 : 0.10);
        applyFlowArrow(sel, isFocus, flowDirection(e.rxn_id, e.met_id), flowMagnitude(e.rxn_id, e.met_id));
      });

    d3.select(zl).select('.g-rail').selectAll<SVGPathElement, RailRouteGeom>('path.rail-line')
      .each(function(d: RailRouteGeom) {
        const rv = railViewById.get(d.line_id);
        const isFocus = active && focusState!.rxnIds.has(d.rxn_id) && focusState!.metIds.has(d.met_id);
        const sel = d3.select(this);
        sel
          .attr('stroke', !active || isFocus ? (rv?.color ?? 'rgba(110,118,129,0.18)') : dimStroke)
          .attr('stroke-width', rv?.width ?? 0.8)
          .attr('opacity', !active || isFocus ? 1 : 0.12);
        applyFlowArrow(sel, isFocus, flowDirection(d.rxn_id, d.met_id, d.stoich), flowMagnitude(d.rxn_id, d.met_id, d.stoich));
      });

    d3.select(zl).select('.g-route').selectAll<SVGPathElement, {rxn_ids?: string[]; pair_id:string; lineIndex:number}>('path.route-line')
      .each(function(d) {
        const sv = view.stations.find(s => s.pair_id === d.pair_id);
        const lv = sv ? sv.lines[d.lineIndex] : undefined;
        const isFocus = active && (d.rxn_ids ?? []).some(rid => focusState!.rxnIds.has(rid));
        d3.select(this)
          .classed('flov-flow-arrow', false)
          .classed('flov-flow-forward', false)
          .classed('flov-flow-reverse', false)
          .style('animation-duration', '')
          .attr('marker-end', null)
          .attr('marker-start', null)
          .attr('stroke', !active || isFocus ? (lv?.color ?? '#6e7681') : dimStroke)
          .attr('stroke-width', lv?.width ?? 1)
          .attr('opacity', !active || isFocus ? 1 : 0.10);
      });

    d3.select(zl).select('.g-stn').selectAll<SVGRectElement, {st: StationGeom; side: 'a' | 'b'}>('rect.stn-pill')
      .each(function(p) {
        const stroke = p.side === 'a' ? p.st.comp_a_color : p.st.comp_b_color;
        const isFocus = active && focusHasStation(p.st);
        d3.select(this)
          .attr('fill', !active || isFocus ? '#ffffff' : dimFill)
          .attr('stroke', !active || isFocus ? stroke : dimStroke)
          .attr('opacity', !active || isFocus ? 1 : 0.16);
      });

    gMetSel?.selectAll<SVGPathElement, MetNode>('path.met-shape')
      .each(function(d) {
        const isFocus = active && focusState!.metIds.has(d.id);
        d3.select(this)
          .attr('fill', !active || isFocus ? (d.display_kind === 'stop' ? '#ffffff' : d.comp_color) : dimFill)
          .attr('fill-opacity', !active || isFocus ? (d.display_kind === 'stop' ? 1.0 : 0.85) : 0.55)
          .attr('stroke', !active || isFocus ? (d.display_kind === 'stop' ? d.comp_color : '#0d1117') : dimStroke)
          .attr('opacity', !active || isFocus ? 1 : 0.18);
      });

    gRxnSel?.each(function(d) {
      const isFocus = active && focusState!.rxnIds.has(d.id);
      d3.select(this)
        .attr('fill', !active || isFocus ? d.fill_color : dimFill)
        .attr('stroke', !active || isFocus ? d.border_color : dimStroke)
        .attr('stroke-width', d.border_width)
        .attr('opacity', !active || isFocus ? d.opacity : 0.16);
    });

    d3.select(zl).select('.g-lbl').selectAll<SVGTextElement, RxnNode>('text')
      .attr('opacity', d => !active || focusState!.rxnIds.has(d.id) ? 1 : 0.10);
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

    // Dark/light mode
    const themeGroup = btnGroup(bar, 'Theme');
    const themeBtn = mkBtn('Dark mode', true);
    const applyTheme = (): void => {
      wrapper.classList.toggle('light-mode', lightMode);
      if (svgEl) svgEl.style.background = lightMode ? '#f6f8fa' : '#0d1117';
      themeBtn.textContent = lightMode ? 'Dark mode' : 'Light mode';
      themeBtn.classList.toggle('active', lightMode);
      applyFocus();
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
    d3.select(zl).select('.g-route').selectAll<SVGPathElement, RouteDatum>('path.route-line')
      .each(function(d: RouteDatum) {
        const sv = view.stations.find(s => s.pair_id === d.pair_id);
        const lv = sv ? sv.lines[d.lineIndex] : undefined;
        if (lv) {
          d3.select(this)
            .attr('d', railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k))
            .attr('stroke', lv.color)
            .attr('stroke-width', lv.width);
        } else {
          d3.select(this)
            .attr('d', railOffsetPath(d.path, d.seg_slots, xS, yS, currentTransform.k))
            .attr('stroke', '#6e7681')
            .attr('stroke-width', 1);
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
    applyFocus();
    refreshStationPanelForView();
  }

  // ── Station panel ──

  function stationClassLabel(klass: string): string {
    return klass.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
  }

  function compactMetList(mets: string[] | undefined): string {
    if (!mets || !mets.length) return '---';
    const text = mets.slice(0, 2).join(', ');
    return mets.length > 2 ? `${text}, ...` : text;
  }

  function activeFluxLabel(rxn: RxnSummary): string {
    const raw = rxn.flux_per_view[currentView] ?? rxn.flux_per_view[0] ?? '';
    const idx = raw.indexOf(':');
    return idx >= 0 ? raw.slice(idx + 1).trim() : raw;
  }

  function stationReactionSearchText(rxn: RxnSummary): string {
    return [
      rxn.id,
      rxn.name,
      ...(rxn.substrates ?? []),
      ...(rxn.products ?? []),
      ...(rxn.flux_per_view ?? []),
    ].join(' ').toLowerCase();
  }

  function buildStationPanel(parent: HTMLElement): void {
    const el = document.createElement('div');
    el.className = 'flov-station-panel';
    el.addEventListener('click', ev => ev.stopPropagation());
    parent.appendChild(el);
    stationPanelEl = el;
  }

  function hideStationPanel(): void {
    if (!stationPanelEl) return;
    stationPanelEl.style.display = 'none';
    stationPanelEl.innerHTML = '';
    activeStationPanel = null;
  }

  function refreshStationPanelForView(): void {
    if (!stationPanelEl || stationPanelEl.style.display === 'none' || !activeStationPanel) {
      return;
    }
    showStationPanel(activeStationPanel, true);
  }

  function showStationPanel(st: StationGeom, preserveState = false): void {
    if (!stationPanelEl) return;
    const scrollTop = preserveState ? stationPanelEl.scrollTop : 0;
    const openGroups = preserveState
      ? new Set(
          Array.from(stationPanelEl.querySelectorAll<HTMLElement>('.flov-station-group.open'))
            .map(group => group.dataset.stationClass ?? '')
            .filter(Boolean)
        )
      : new Set<string>();
    activeStationPanel = st;
    stationPanelEl.innerHTML = '';
    stationPanelEl.style.borderColor = st.comp_a_color;

    const head = document.createElement('div');
    head.className = 'flov-station-head';
    const title = document.createElement('div');
    title.className = 'flov-station-title';
    const titleText = document.createElement('span');
    titleText.textContent = `${st.comp_a_label} ↔ ${st.comp_b_label}`;
    const count = document.createElement('span');
    count.className = 'flov-station-badge';
    count.textContent = `${st.n_rxns} rxn${st.n_rxns === 1 ? '' : 's'}`;
    const close = document.createElement('button');
    close.className = 'flov-station-close';
    close.type = 'button';
    close.setAttribute('aria-label', 'Close station flows');
    close.textContent = '×';
    close.onclick = ev => {
      ev.stopPropagation();
      if (focusState?.type === 'station' && focusState.id === st.pair_id) {
        clearFocus();
      } else {
        hideStationPanel();
      }
    };
    title.appendChild(titleText);
    title.appendChild(count);
    title.appendChild(close);
    const sub = document.createElement('div');
    sub.className = 'flov-station-sub';
    sub.textContent = `Reaction flows · ${fig?.views[currentView]?.label ?? ''}`;
    head.appendChild(title);
    head.appendChild(sub);

    const body = document.createElement('div');
    body.className = 'flov-station-body';
    st.lines.forEach(ln => {
      const group = document.createElement('div');
      group.className = 'flov-station-group';
      group.dataset.stationClass = ln.klass;
      const toggle = document.createElement('button');
      toggle.className = 'flov-station-toggle';
      toggle.type = 'button';
      const caret = document.createElement('span');
      caret.className = 'flov-station-caret';
      caret.textContent = '▸';
      const swatch = document.createElement('span');
      swatch.className = 'flov-station-swatch';
      swatch.style.background = KLASS_COLORS[ln.klass] ?? '#6e7681';
      const label = document.createElement('span');
      label.textContent = stationClassLabel(ln.klass);
      const rxnCount = document.createElement('span');
      rxnCount.className = 'flov-station-count';
      rxnCount.textContent = `${ln.rxn_summaries.length}`;
      toggle.appendChild(caret);
      toggle.appendChild(swatch);
      toggle.appendChild(label);
      toggle.appendChild(rxnCount);

      const rxns = document.createElement('div');
      rxns.className = 'flov-station-rxns';
      appendStationRailMaps(rxns, ln, st);

      toggle.onclick = ev => {
        ev.stopPropagation();
        const open = group.classList.toggle('open');
        caret.textContent = open ? '▾' : '▸';
      };
      group.appendChild(toggle);
      group.appendChild(rxns);
      body.appendChild(group);
    });

    stationPanelEl.appendChild(head);
    stationPanelEl.appendChild(body);
    stationPanelEl.style.display = 'block';
    if (preserveState) {
      stationPanelEl.querySelectorAll<HTMLElement>('.flov-station-group').forEach(group => {
        if (!openGroups.has(group.dataset.stationClass ?? '')) return;
        group.classList.add('open');
        const caret = group.querySelector<HTMLElement>('.flov-station-caret');
        if (caret) caret.textContent = '▾';
      });
      stationPanelEl.scrollTop = scrollTop;
    } else {
      stationPanelEl.scrollTop = 0;
    }
    applyStationPanelSearch(currentStationSearch);
  }

  function applyStationPanelSearch(query: string, scrollToRxnId?: string): void {
    if (!stationPanelEl || stationPanelEl.style.display === 'none') return;
    const lq = query.trim().toLowerCase();
    const shouldFilter = lq.length > 0;

    stationPanelEl.querySelectorAll<HTMLElement>('.flov-station-group').forEach(group => {
      const branches = Array.from(group.querySelectorAll<SVGGElement>('.flov-station-branch'));
      let groupHasMatch = false;
      branches.forEach(branch => {
        const isTarget = scrollToRxnId !== undefined && branch.dataset.rxnId === scrollToRxnId;
        const isMatch = !shouldFilter || (branch.dataset.searchText ?? '').includes(lq);
        groupHasMatch ||= isMatch || isTarget;
        branch.classList.toggle('search-hidden', shouldFilter && !isMatch && !isTarget);
        branch.classList.toggle('search-hit', (shouldFilter && isMatch) || isTarget);
      });
      group.querySelectorAll<SVGSVGElement>('.flov-station-map').forEach(svg => {
        const hasVisibleBranch = Array.from(
          svg.querySelectorAll<SVGGElement>('.flov-station-branch')
        ).some(branch => !branch.classList.contains('search-hidden'));
        svg.classList.toggle('search-hidden', shouldFilter && !hasVisibleBranch);
      });

      group.classList.toggle('search-hidden', shouldFilter && !groupHasMatch);
      const caret = group.querySelector<HTMLElement>('.flov-station-caret');
      if (shouldFilter && groupHasMatch) {
        group.classList.add('open');
        if (caret) caret.textContent = '▾';
      } else if (!shouldFilter) {
        branches.forEach(branch => branch.classList.remove('search-hit', 'search-hidden'));
        group.querySelectorAll<SVGSVGElement>('.flov-station-map').forEach(svg => {
          svg.classList.remove('search-hidden');
        });
        group.classList.remove('search-hidden');
      }
    });

    if (scrollToRxnId) {
      const target = Array.from(
        stationPanelEl.querySelectorAll<SVGGElement>('.flov-station-branch')
      ).find(branch => branch.dataset.rxnId === scrollToRxnId);
      target?.scrollIntoView({ block: 'center', inline: 'nearest' });
    }
  }

  function appendStationRailMaps(
    parent: HTMLElement,
    line: StationLineGeom,
    station: StationGeom
  ): void {
    const svgNS = 'http://www.w3.org/2000/svg';
    const color = KLASS_COLORS[line.klass] ?? '#6e7681';
    const chunkSize = 5;
    for (let start = 0; start < line.rxn_summaries.length; start += chunkSize) {
      const chunk = line.rxn_summaries.slice(start, start + chunkSize);
      const width = 650;
      const height = 225;
      const railY = 160;
      const x0 = 112;
      const splitStartX = 205;
      const splitGap = 62;
      const x1 = width - 126;
      const laneGap = 9;
      const markerPad = 4;
      const laneTop = railY - ((chunk.length - 1) / 2) * laneGap;
      const laneBottom = railY + ((chunk.length - 1) / 2) * laneGap;
      const bundleTop = laneTop - markerPad;
      const bundleBottom = laneBottom + markerPad;
      const svg = document.createElementNS(svgNS, 'svg');
      svg.setAttribute('class', 'flov-station-map');
      svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
      svg.setAttribute('role', 'img');
      svg.setAttribute(
        'aria-label',
        `${stationClassLabel(line.klass)} reaction rail ${start / chunkSize + 1}`
      );

      const endLaneY = railY + (chunk.length - 1 - (chunk.length - 1) / 2) * laneGap;
      [
        [x0, station.comp_a_color, bundleTop, bundleBottom],
        [x1, station.comp_b_color, endLaneY - markerPad, endLaneY + markerPad],
      ].forEach(([x, stroke, yTop, yBottom]) => {
        const pill = document.createElementNS(svgNS, 'rect');
        pill.setAttribute('x', String(Number(x) - 4));
        pill.setAttribute('y', String(yTop));
        pill.setAttribute('width', '8');
        pill.setAttribute('height', String(Number(yBottom) - Number(yTop)));
        pill.setAttribute('rx', '4');
        pill.setAttribute('fill', '#fff');
        pill.setAttribute('stroke', String(stroke));
        pill.setAttribute('stroke-width', '2');
        svg.appendChild(pill);
      });

      chunk.forEach((rxn, i) => {
        const isThroughLine = i === chunk.length - 1;
        const laneY = railY + (i - (chunk.length - 1) / 2) * laneGap;
        const splitX = splitStartX + i * splitGap;
        const branchEndX = isThroughLine
          ? x1 - 34
          : Math.min(splitX + 76, x1 - 42);
        const branchEndY = isThroughLine ? laneY : 34 + i * 31;
        const labelX = isThroughLine
          ? (splitStartX + Math.max(0, chunk.length - 2) * splitGap + x1) / 2 + 42
          : (branchEndX + x1) / 2;
        const labelY = branchEndY - (isThroughLine ? 28 : 14);
        const fluxX = isThroughLine
          ? (splitStartX + Math.max(0, chunk.length - 2) * splitGap + x1) / 2
          : ((i === 0 ? x0 : splitStartX + (i - 1) * splitGap) + splitX) / 2;
        const fluxY = laneY - 5;
        const labelAnchor = branchEndX > width * 0.76 ? 'end' : 'middle';
        const branchColor = STATION_BRANCH_COLORS[i % STATION_BRANCH_COLORS.length];
        const substrate = compactMetList(rxn.substrates);
        const product = compactMetList(rxn.products);

        const g = document.createElementNS(svgNS, 'g');
        g.setAttribute('class', 'flov-station-branch');
        g.dataset.rxnId = rxn.id;
        g.dataset.searchText = stationReactionSearchText(rxn);
        g.addEventListener('click', ev => setStationReactionFocus(rxn.id, ev as MouseEvent));

        const subLabel = document.createElementNS(svgNS, 'text');
        subLabel.setAttribute('class', 'flov-station-map-met');
        subLabel.setAttribute('x', String(x0 - 10));
        subLabel.setAttribute('y', (laneY + 3).toFixed(1));
        subLabel.setAttribute('text-anchor', 'end');
        subLabel.textContent = substrate.length > 20 ? `${substrate.slice(0, 19)}...` : substrate;
        g.appendChild(subLabel);

        const branch = document.createElementNS(svgNS, 'path');
        branch.setAttribute(
          'd',
          isThroughLine
            ? `M${x0} ${laneY.toFixed(1)} L${x1} ${laneY.toFixed(1)}`
            : `M${x0} ${laneY.toFixed(1)} L${splitX} ${laneY.toFixed(1)} L${branchEndX.toFixed(1)} ${branchEndY.toFixed(1)} L${x1} ${branchEndY.toFixed(1)}`
        );
        branch.setAttribute('class', 'flov-station-map-line');
        branch.setAttribute('stroke', branchColor);
        branch.setAttribute('stroke-width', '1.9');
        branch.setAttribute('fill', 'none');
        g.appendChild(branch);

        const productLabel = document.createElementNS(svgNS, 'text');
        productLabel.setAttribute('class', 'flov-station-map-met');
        productLabel.setAttribute('x', String(x1 + 10));
        productLabel.setAttribute('y', (branchEndY + 3).toFixed(1));
        productLabel.setAttribute('text-anchor', 'start');
        productLabel.textContent = product.length > 20 ? `${product.slice(0, 19)}...` : product;
        g.appendChild(productLabel);

        if (!isThroughLine) {
          const dot = document.createElementNS(svgNS, 'circle');
          dot.setAttribute('cx', String(splitX));
          dot.setAttribute('cy', laneY.toFixed(1));
          dot.setAttribute('r', '2.7');
          dot.setAttribute('fill', '#fff');
          dot.setAttribute('stroke', branchColor);
          dot.setAttribute('stroke-width', '1.5');
          g.appendChild(dot);
        }

        const label = document.createElementNS(svgNS, 'text');
        label.setAttribute('class', 'flov-station-map-label');
        label.setAttribute('x', labelX.toFixed(1));
        label.setAttribute('y', String(labelY));
        label.setAttribute('text-anchor', labelAnchor);
        label.textContent = rxn.id.length > 18 ? `${rxn.id.slice(0, 17)}...` : rxn.id;
        g.appendChild(label);

        const flux = document.createElementNS(svgNS, 'text');
        flux.setAttribute('class', 'flov-station-map-flux');
        flux.setAttribute('x', fluxX.toFixed(1));
        flux.setAttribute('y', String(fluxY));
        flux.setAttribute('text-anchor', 'middle');
        flux.textContent = activeFluxLabel(rxn);
        g.appendChild(flux);

        const title = document.createElementNS(svgNS, 'title');
        title.textContent = `${rxn.id}${rxn.name && rxn.name !== rxn.id ? ` · ${rxn.name}` : ''}\nSub: ${compactMetList(rxn.substrates)}\nPrd: ${compactMetList(rxn.products)}\n${rxn.flux_per_view.join('\n')}`;
        g.appendChild(title);
        svg.appendChild(g);
      });
      chunk.slice(0, -1).forEach((_, i) => {
        const splitX = splitStartX + i * splitGap;
        const remainingLaneYs = chunk
          .slice(i)
          .map((__, ri) => railY + (i + ri - (chunk.length - 1) / 2) * laneGap);
        const markerTop = Math.min(...remainingLaneYs) - markerPad;
        const markerBottom = Math.max(...remainingLaneYs) + markerPad;
        const marker = document.createElementNS(svgNS, 'rect');
        marker.setAttribute('x', String(splitX - 3.5));
        marker.setAttribute('y', String(markerTop));
        marker.setAttribute('width', '7');
        marker.setAttribute('height', String(markerBottom - markerTop));
        marker.setAttribute('rx', '3.5');
        marker.setAttribute('fill', '#fff');
        marker.setAttribute('stroke', color);
        marker.setAttribute('stroke-width', '1.5');
        svg.appendChild(marker);
      });
      parent.appendChild(svg);
    }
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

  // ── Legend ──

  function buildLegend(parent: HTMLElement, data: D3FigureData, svgH: number): void {
    const SVG_NS = 'http://www.w3.org/2000/svg';
    const gradId = `flov-legend-grad-${Math.random().toString(36).slice(2)}`;

    const panel = document.createElement('div');
    panel.className = 'flov-legend';
    panel.style.maxHeight = Math.max(svgH - 70, 280) + 'px';
    parent.appendChild(panel);

    const head = document.createElement('div');
    head.className = 'flov-legend-head';
    const headTitle = document.createElement('span');
    headTitle.className = 'flov-legend-title';
    headTitle.textContent = 'Legend';
    const headCaret = document.createElement('span');
    headCaret.className = 'flov-legend-caret';
    headCaret.textContent = '▾';
    head.appendChild(headTitle);
    head.appendChild(headCaret);
    panel.appendChild(head);

    const body = document.createElement('div');
    body.className = 'flov-legend-body';
    panel.appendChild(body);

    head.onclick = () => {
      const open = !panel.classList.contains('collapsed');
      panel.classList.toggle('collapsed', open);
      headCaret.textContent = open ? '▸' : '▾';
    };

    function section(title: string): HTMLDivElement {
      const sec = document.createElement('div');
      sec.className = 'flov-legend-section';
      const t = document.createElement('div');
      t.className = 'flov-legend-section-title';
      t.textContent = title;
      sec.appendChild(t);
      body.appendChild(sec);
      return sec;
    }

    function row(sec: HTMLDivElement, swatch: SVGSVGElement, label: string, sub?: string): void {
      const r = document.createElement('div');
      r.className = 'flov-legend-row';
      const sw = document.createElement('div');
      sw.className = 'flov-legend-swatch';
      sw.appendChild(swatch);
      const text = document.createElement('div');
      text.className = 'flov-legend-text';
      const main = document.createElement('div');
      main.className = 'flov-legend-label';
      main.textContent = label;
      text.appendChild(main);
      if (sub) {
        const s = document.createElement('div');
        s.className = 'flov-legend-sub';
        s.textContent = sub;
        text.appendChild(s);
      }
      r.appendChild(sw);
      r.appendChild(text);
      sec.appendChild(r);
    }

    function mkSvg(w: number, h: number): SVGSVGElement {
      const s = document.createElementNS(SVG_NS, 'svg') as SVGSVGElement;
      s.setAttribute('width', String(w));
      s.setAttribute('height', String(h));
      s.setAttribute('viewBox', `0 0 ${w} ${h}`);
      return s;
    }
    function mkEl(tag: string, attrs: Record<string, string | number>): SVGElement {
      const e = document.createElementNS(SVG_NS, tag);
      for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, String(v));
      return e;
    }

    // ── Flux scale (horizontal gradient bar) ──
    {
      const sec = section('Flux (mmol/gDW/h)');
      const cbW = 198, cbH = 38;
      const svg = mkSvg(cbW, cbH);
      const defs = mkEl('defs', {});
      const grad = mkEl('linearGradient', { id: gradId, x1: '0', y1: '0', x2: '1', y2: '0' });
      // left → right: dark blue → light blue → grey → light red → dark red
      const stops: [string, string][] = [
        ['0%',   '#1565c0'],
        ['25%',  '#90caf9'],
        ['50%',  '#6e7681'],
        ['75%',  '#ef9a9a'],
        ['100%', '#c62828'],
      ];
      stops.forEach(([off, col]) => {
        grad.appendChild(mkEl('stop', { offset: off, 'stop-color': col }));
      });
      defs.appendChild(grad);
      svg.appendChild(defs);
      svg.appendChild(mkEl('rect', { x: 2, y: 4, width: cbW - 4, height: 12, rx: 3, fill: `url(#${gradId})` }));

      const absMax = data.meta.abs_max_flux;
      const fmt = d3.format('+.2f');
      const labels: [number, string][] = [
        [2,        fmt(-absMax)],
        [cbW / 2,  '0'],
        [cbW - 2,  fmt(absMax)],
      ];
      labels.forEach(([x, txt], i) => {
        const t = mkEl('text', {
          x, y: 28,
          'text-anchor': i === 0 ? 'start' : i === labels.length - 1 ? 'end' : 'middle',
          'font-size': 9,
          'font-family': "'Inter',Arial,sans-serif",
          fill: '#8b949e',
        });
        t.textContent = txt;
        svg.appendChild(t);
      });

      const wrap = document.createElement('div');
      wrap.className = 'flov-legend-fluxbar';
      wrap.appendChild(svg);
      sec.appendChild(wrap);

      const note = document.createElement('div');
      note.className = 'flov-legend-note';
      note.textContent = 'Grey = inactive / no data';
      sec.appendChild(note);
    }

    // ── Reactions ──
    {
      const sec = section('Reactions');
      // Circle with flux fill + compartment border (the standard node)
      const sv1 = mkSvg(28, 20);
      sv1.appendChild(mkEl('circle', { cx: 14, cy: 10, r: 7, fill: '#ef9a9a', stroke: '#58a6ff', 'stroke-width': 2 }));
      row(sec, sv1, 'Reaction node', 'Fill = flux  ·  Border = compartment');

      // Faded (no flux data)
      const sv2 = mkSvg(28, 20);
      sv2.appendChild(mkEl('circle', { cx: 14, cy: 10, r: 7, fill: '#6e7681', stroke: '#58a6ff', 'stroke-width': 2, opacity: 0.22 }));
      row(sec, sv2, 'No flux data', 'Faded grey fill');
    }

    // ── Metabolites ──
    {
      const sec = section('Metabolites');
      // Stop: white circle with comp border (routed)
      const sv1 = mkSvg(28, 20);
      sv1.appendChild(mkEl('circle', { cx: 14, cy: 10, r: 6.5, fill: '#ffffff', stroke: '#58a6ff', 'stroke-width': 1.6 }));
      row(sec, sv1, 'Routed metabolite', 'White stop on a rail line');

      // Detached: filled diamond
      const sv2 = mkSvg(28, 20);
      sv2.appendChild(mkEl('path', {
        d: 'M14,3 L21,10 L14,17 L7,10 Z',
        fill: '#58a6ff', 'fill-opacity': 0.85, stroke: '#0d1117', 'stroke-width': 0.7,
      }));
      row(sec, sv2, 'Detached metabolite', 'Off-rail, no station');
    }

    // ── Compartments (dynamic from data) ──
    if (data.compartments && data.compartments.length) {
      const sec = section('Compartments');
      for (const c of data.compartments) {
        const sv = mkSvg(28, 20);
        sv.appendChild(mkEl('rect', {
          x: 3, y: 4, width: 22, height: 12, rx: 3,
          fill: c.fill, stroke: c.color, 'stroke-width': 1.4,
          'stroke-dasharray': '4,2',
        }));
        row(sec, sv, c.label || c.key);
      }
    }

    // ── Stations ──
    {
      const sec = section('Stations');
      // Pill: rotated rounded rectangle, two compartment-colored ends
      const sv = mkSvg(28, 28);
      const g = mkEl('g', { transform: 'translate(14,14)' });
      // Two stacked pills representing the two compartment sides (a top, b bottom)
      g.appendChild(mkEl('rect', {
        x: -4, y: -10, width: 8, height: 10, rx: 4, ry: 4,
        fill: '#ffffff', stroke: '#3fb950', 'stroke-width': 1.8,
      }));
      g.appendChild(mkEl('rect', {
        x: -4, y: 0, width: 8, height: 10, rx: 4, ry: 4,
        fill: '#ffffff', stroke: '#58a6ff', 'stroke-width': 1.8,
      }));
      sv.appendChild(g);
      row(sec, sv, 'Station hub', 'Hub between two compartments');
    }

    // ── Station classes (route line colours) ──
    {
      const sec = section('Station line classes');
      const entries: [string, string][] = [
        ['amino_acid', 'Amino acids'],
        ['sugar',      'Sugars'],
        ['cofactor',   'Cofactors'],
        ['inorganic',  'Inorganics'],
        ['other',      'Other'],
      ];
      for (const [k, lbl] of entries) {
        const sv = mkSvg(28, 14);
        sv.appendChild(mkEl('line', {
          x1: 2, y1: 7, x2: 26, y2: 7,
          stroke: KLASS_COLORS[k] ?? '#ec4899',
          'stroke-width': 3.4,
          'stroke-linecap': 'round',
        }));
        row(sec, sv, lbl);
      }
      // Inactive
      const svI = mkSvg(28, 14);
      svI.appendChild(mkEl('line', {
        x1: 2, y1: 7, x2: 26, y2: 7,
        stroke: '#6e7681', 'stroke-width': 2.2, 'stroke-linecap': 'round',
      }));
      row(sec, svI, 'Inactive class', 'No flux on this line');
    }

    // ── Line types ──
    {
      const sec = section('Lines');
      // Solid rail line (rxn ↔ met, active)
      const sv1 = mkSvg(28, 14);
      sv1.appendChild(mkEl('line', {
        x1: 2, y1: 7, x2: 26, y2: 7, stroke: '#58a6ff',
        'stroke-width': 2.6, 'stroke-linecap': 'round',
      }));
      row(sec, sv1, 'Active line', 'Solid: flux on this edge');

      // Striped forward flow (rxn → met)
      const sv2 = mkSvg(28, 14);
      const ln2 = mkEl('line', {
        x1: 2, y1: 7, x2: 26, y2: 7, stroke: '#58a6ff',
        'stroke-width': 2.6, 'stroke-linecap': 'round',
      });
      ln2.setAttribute('class', 'flov-flow-arrow flov-flow-forward');
      sv2.appendChild(ln2);
      row(sec, sv2, 'Flow: rxn → met', 'Forward direction (striped)');

      // Striped reverse flow (met → rxn)
      const sv3 = mkSvg(28, 14);
      const ln3 = mkEl('line', {
        x1: 2, y1: 7, x2: 26, y2: 7, stroke: '#58a6ff',
        'stroke-width': 2.6, 'stroke-linecap': 'round',
      });
      ln3.setAttribute('class', 'flov-flow-arrow flov-flow-reverse');
      sv3.appendChild(ln3);
      row(sec, sv3, 'Flow: met → rxn', 'Reverse direction (striped)');

      // Inactive grey line
      const sv4 = mkSvg(28, 14);
      sv4.appendChild(mkEl('line', {
        x1: 2, y1: 7, x2: 26, y2: 7, stroke: '#6e7681',
        'stroke-width': 2.2, 'stroke-linecap': 'round', opacity: 0.6,
      }));
      row(sec, sv4, 'Inactive', 'Grey: no flux on this edge');
    }
  }

  // ── Search ──

  function buildSearch(parent: HTMLElement, data: D3FigureData): void {
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

    function selectResult(type: 'rxn' | 'met', idx: number): void {
      inp.value = type === 'rxn'
        ? (ia.rxn_names[idx] || ia.rxn_ids[idx])
        : (ia.met_names[idx] || ia.met_ids[idx]);
      res.style.display = 'none';
      if (type === 'rxn') {
        focusStationReactionById(ia.rxn_ids[idx]);
      } else {
        focusMetaboliteById(ia.met_ids[idx]);
      }
    }

    type StationHit = {
      type: 'station-rxn';
      station: StationGeom;
      line: StationLineGeom;
      rxn: RxnSummary;
    };

    function selectStationResult(hit: StationHit): void {
      showStationPanel(hit.station);
      focusStationReactionById(hit.rxn.id);
      inp.value = hit.rxn.name && hit.rxn.name !== hit.rxn.id ? hit.rxn.name : hit.rxn.id;
      currentStationSearch = inp.value;
      applyStationPanelSearch(currentStationSearch, hit.rxn.id);
      res.style.display = 'none';
    }

    function runSearch(q: string): void {
      res.innerHTML = '';
      currentStationSearch = q;
      applyStationPanelSearch(q);
      if (!q.trim()) { res.style.display = 'none'; return; }
      const lq = q.toLowerCase();
      const hits: ({ type: 'rxn' | 'met'; idx: number } | StationHit)[] = [];

      for (let i = 0; i < ia.rxn_ids.length; i++) {
        if (rxnIdsLc[i].includes(lq) || rxnNamesLc[i].includes(lq)) {
          if (hits.length < 20) hits.push({ type: 'rxn', idx: i });
        }
      }
      for (let i = 0; i < ia.met_ids.length; i++) {
        if (metIdsLc[i].includes(lq) || metNamesLc[i].includes(lq)) {
          if (hits.length < 20) hits.push({ type: 'met', idx: i });
        }
      }
      for (const station of data.stations) {
        for (const line of station.lines) {
          for (const rxn of line.rxn_summaries) {
            if (!stationReactionSearchText(rxn).includes(lq)) continue;
            if (hits.length < 20) hits.push({ type: 'station-rxn', station, line, rxn });
          }
        }
      }

      if (!hits.length) {
        const d = document.createElement('div');
        d.className = 'flov-search-result';
        d.style.color = '#7d8590';
        d.textContent = 'No matches';
        res.appendChild(d);
        res.style.display = 'block';
        return;
      }

      for (const h of hits) {
        const el = document.createElement('div');
        el.className = 'flov-search-result';
        if (h.type === 'station-rxn') {
          const id = h.rxn.id;
          const name = h.rxn.name;
          el.title = `${id}${name && name !== id ? ' — ' + name : ''}`;
          el.innerHTML =
            `<b>${id}</b>` +
            (name && name !== id ? ` <span style="color:#7d8590">— ${name}</span>` : '') +
            ` <span style="color:#444c56;font-size:10px">[station ${stationClassLabel(h.line.klass)}]</span>`;
          el.onclick = () => selectStationResult(h);
        } else {
          const id   = h.type === 'rxn' ? ia.rxn_ids[h.idx]   : ia.met_ids[h.idx];
          const name = h.type === 'rxn' ? ia.rxn_names[h.idx]  : ia.met_names[h.idx];
          el.title = id + (name && name !== id ? ' — ' + name : '');
          el.innerHTML =
            `<b>${id}</b>` +
            (name && name !== id ? ` <span style="color:#7d8590">— ${name}</span>` : '') +
            ` <span style="color:#444c56;font-size:10px">[${h.type}]</span>`;
          el.onclick = () => selectResult(h.type, h.idx);
        }
        res.appendChild(el);
      }
      res.style.display = 'block';
    }

    inp.addEventListener('input', function (this: HTMLInputElement) { runSearch(this.value); });
    inp.addEventListener('keydown', (e: KeyboardEvent) => {
      if (e.key === 'Escape') {
        inp.value = '';
        currentStationSearch = '';
        applyStationPanelSearch('');
        if (focusState) {
          focusState = null;
          applyFocus();
        }
        res.style.display = 'none';
      }
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
