import * as d3 from 'd3';
import type { HoverData, MetNode } from './types';

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
.flov-station-panel{position:absolute;top:52px;left:10px;z-index:95;width:min(720px,calc(100% - 20px));max-height:min(520px,calc(100% - 70px));overflow-y:auto;background:#161b22;border:1px solid #30363d;border-radius:8px;box-shadow:0 10px 28px rgba(0,0,0,0.62);color:#e6edf3;font:12px/1.45 'Inter',Arial,sans-serif;display:none;}
.flov-station-head{position:sticky;top:0;background:#161b22;border-bottom:1px solid #30363d;padding:10px 12px;z-index:1;}
.flov-station-title{display:flex;align-items:center;gap:8px;font-weight:700;color:#fff;}
.flov-station-badge{background:#21262d;border-radius:3px;padding:1px 5px;font-size:10px;color:#8b949e;letter-spacing:.04em;}
.flov-station-close{margin-left:auto;background:transparent;border:0;color:#8b949e;font:18px/1 Arial,sans-serif;cursor:pointer;padding:0 2px;}
.flov-station-close:hover{color:#e6edf3;}
.flov-station-sub{margin-top:3px;color:#8b949e;font-size:11px;}
.flov-station-body{padding:7px 8px 10px;}
.flov-station-group{border-top:1px solid #21262d;}
.flov-station-group:first-child{border-top:0;}
.flov-station-toggle{display:flex;align-items:center;gap:7px;width:100%;padding:8px 5px;background:transparent;border:0;color:#e6edf3;text-align:left;font:12px/1.2 'Inter',Arial,sans-serif;cursor:pointer;}
.flov-station-toggle:hover{background:#21262d;}
.flov-station-caret{width:12px;color:#8b949e;font-size:11px;}
.flov-station-swatch{width:9px;height:9px;border-radius:2px;flex:0 0 auto;}
.flov-station-count{margin-left:auto;color:#7d8590;font-size:10px;}
.flov-station-rxns{display:none;padding:0 5px 8px 31px;}
.flov-station-group.open .flov-station-rxns{display:block;}
.flov-station-map{width:100%;height:auto;margin:5px 0 10px;overflow:visible;}
.flov-station-branch{cursor:pointer;}
.flov-station-branch:hover .flov-station-map-line{stroke-width:2.4;}
.flov-station-branch.search-hit .flov-station-map-line{stroke-width:3.1;}
.flov-station-branch.search-hidden{display:none;}
.flov-station-map.search-hidden{display:none;}
.flov-station-group.search-hidden{display:none;}
.flov-station-map-line{stroke-linecap:round;stroke-linejoin:round;}
.flov-station-map-label{font:12px ui-monospace,SFMono-Regular,Menlo,monospace;fill:#8b949e;}
.flov-station-map-flux{font:10px ui-monospace,SFMono-Regular,Menlo,monospace;fill:#7d8590;}
.flov-station-map-met{font:9px 'Inter',Arial,sans-serif;fill:#8b949e;}
.flov-legend{position:absolute;top:52px;right:8px;z-index:90;width:230px;overflow-y:auto;background:#161b22;border:1px solid #30363d;border-radius:8px;box-shadow:0 4px 14px rgba(0,0,0,0.55);color:#e6edf3;font:11px/1.45 'Inter',Arial,sans-serif;}
.flov-legend-head{display:flex;align-items:center;gap:6px;padding:7px 10px;border-bottom:1px solid #30363d;cursor:pointer;user-select:none;background:#1c2128;border-top-left-radius:8px;border-top-right-radius:8px;}
.flov-legend-title{font-weight:700;color:#fff;flex:1;font-size:12px;letter-spacing:.02em;}
.flov-legend-caret{color:#8b949e;font-size:11px;}
.flov-legend-body{padding:4px 10px 10px;}
.flov-legend.collapsed .flov-legend-body{display:none;}
.flov-legend.collapsed{border-bottom-left-radius:8px;border-bottom-right-radius:8px;}
.flov-legend.collapsed .flov-legend-head{border-bottom:0;}
.flov-legend-section{margin:8px 0 4px;}
.flov-legend-section + .flov-legend-section{border-top:1px solid #21262d;padding-top:6px;}
.flov-legend-section-title{font-size:9px;letter-spacing:.08em;text-transform:uppercase;color:#8b949e;margin-bottom:4px;font-weight:600;}
.flov-legend-row{display:flex;align-items:center;gap:8px;margin:3px 0;}
.flov-legend-swatch{flex:0 0 30px;display:flex;justify-content:center;align-items:center;}
.flov-legend-text{flex:1;min-width:0;}
.flov-legend-label{color:#e6edf3;font-size:11px;line-height:1.25;}
.flov-legend-sub{color:#7d8590;font-size:9.5px;line-height:1.2;margin-top:1px;}
.flov-legend-fluxbar{margin:2px 0 0;}
.flov-legend-note{color:#7d8590;font-size:9.5px;margin-top:2px;font-style:italic;}
.flov-wrapper.light-mode .flov-legend{background:#ffffff;border-color:#d0d7de;color:#24292f;box-shadow:0 8px 24px rgba(31,35,40,0.12);}
.flov-wrapper.light-mode .flov-legend-head{background:#f6f8fa;border-bottom-color:#d0d7de;}
.flov-wrapper.light-mode .flov-legend-title{color:#24292f;}
.flov-wrapper.light-mode .flov-legend-caret{color:#57606a;}
.flov-wrapper.light-mode .flov-legend-section + .flov-legend-section{border-top-color:#d8dee4;}
.flov-wrapper.light-mode .flov-legend-section-title,.flov-wrapper.light-mode .flov-legend-sub,.flov-wrapper.light-mode .flov-legend-note{color:#57606a;}
.flov-wrapper.light-mode .flov-legend-label{color:#24292f;}
@keyframes flov-flow-forward{from{stroke-dashoffset:0}to{stroke-dashoffset:-28}}
@keyframes flov-flow-reverse{from{stroke-dashoffset:0}to{stroke-dashoffset:28}}
.flov-flow-arrow{stroke-dasharray:9 5;}
.flov-flow-forward{animation:flov-flow-forward 1s linear infinite;}
.flov-flow-reverse{animation:flov-flow-reverse 1s linear infinite;}
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
.flov-wrapper.light-mode .flov-station-panel{background:#ffffff;color:#24292f;border-color:#d0d7de;box-shadow:0 10px 28px rgba(31,35,40,0.16);}
.flov-wrapper.light-mode .flov-station-head{background:#ffffff;border-bottom-color:#d0d7de;}
.flov-wrapper.light-mode .flov-station-title{color:#24292f;}
.flov-wrapper.light-mode .flov-station-badge{background:#f3f4f6;color:#57606a;}
.flov-wrapper.light-mode .flov-station-close{color:#57606a;}
.flov-wrapper.light-mode .flov-station-close:hover{color:#24292f;}
.flov-wrapper.light-mode .flov-station-sub,.flov-wrapper.light-mode .flov-station-count{color:#57606a;}
.flov-wrapper.light-mode .flov-station-group{border-top-color:#d8dee4;}
.flov-wrapper.light-mode .flov-station-toggle{color:#24292f;}
.flov-wrapper.light-mode .flov-station-toggle:hover{background:#eef2f7;}
.flov-wrapper.light-mode .flov-station-map-label,.flov-wrapper.light-mode .flov-station-map-flux,.flov-wrapper.light-mode .flov-station-map-met{fill:#57606a;}
`;

let _cssInjected = false;
export function injectCSS(): void {
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

export function symPath(symbol: string, r: number): string {
  const area = symArea(r);
  switch (symbol) {
    case 'triangle-up':   return d3.symbol(d3.symbolTriangle, area)(null) ?? '';
    case 'triangle-down': return d3.symbol(d3.symbolTriangle, area)(null) ?? '';
    case 'hexagon':       return d3.symbol(symbolHexagon,      area)(null) ?? '';
    case 'cross':         return d3.symbol(d3.symbolCross,    area * 0.7)(null) ?? '';
    default:              return d3.symbol(d3.symbolCircle,   area)(null) ?? '';
  }
}

export function metPath(m: MetNode): string {
  if (m.display_kind === 'stop') {
    return d3.symbol(d3.symbolCircle, Math.PI * 20)(null) ?? '';
  }
  return d3.symbol(d3.symbolDiamond, Math.PI * 16)(null) ?? '';
}

// ── Path utilities ────────────────────────────────────────────────────────────

export function hullPath(
  verts: [number, number][],
  xS: d3.ScaleLinear<number, number>, yS: d3.ScaleLinear<number, number>
): string {
  const gen = d3.line<[number, number]>()
    .x(d => xS(d[0])).y(d => yS(d[1]))
    .curve(d3.curveCatmullRomClosed.alpha(0.5));
  return gen(verts) ?? '';
}

// Rail path with pixel-constant lane offset.
// Each segment is shifted right-perpendicularly by seg_slots[k] * LANE_W_PX
// divided by the current zoom scale k_zoom, so the screen-pixel gap stays
// constant as the user zooms. Miter intersections keep every segment truly
// parallel; sharp corners fall back to a bevel midpoint.
const LANE_W_PX = 3.5;
const MAX_RENDER_LANE_SLOT = 6.0;

function renderLaneSlot(slot: number): number {
  if (!Number.isFinite(slot)) return 0.0;
  return MAX_RENDER_LANE_SLOT * Math.tanh(slot / MAX_RENDER_LANE_SLOT);
}

export function railOffsetPath(
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
  const off = seg_slots.map(s => (renderLaneSlot(s) * LANE_W_PX) / k_zoom);
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

export function rxnTooltipHTML(h: HoverData): string {
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

export function metTooltipHTML(m: MetNode): string {
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
// flux_network.py. Used for the swatch shown next to each station category.
export const KLASS_COLORS: Record<string, string> = {
  amino_acid: '#22c55e',
  sugar:      '#f59e0b',
  cofactor:   '#a855f7',
  inorganic:  '#06b6d4',
  other:      '#ec4899',
};
export const STATION_BRANCH_COLORS = [
  '#22c55e', '#f59e0b', '#38bdf8', '#ec4899', '#a855f7',
];

