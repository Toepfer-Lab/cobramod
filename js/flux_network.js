function H({ model: S, el: I }) {
  const k = document.createElement("div");
  I.appendChild(k);
  let C = !1, i = null, t = null, c = 0, w = !0, L = null, _ = null, v = null, M = null;
  function N() {
    return window.Plotly ? Promise.resolve() : new Promise((a, d) => {
      const s = document.createElement("script");
      s.src = "https://cdn.plot.ly/plotly-2.35.2.min.js", s.onload = () => a(), s.onerror = () => d(new Error("Failed to load Plotly.js")), document.head.appendChild(s);
    });
  }
  function T(a, d, s) {
    const u = t.met_bounds[a] || [-1e9, 1e9, -1e9, 1e9], [r, x, y, f] = u, p = 0.2;
    let g = r + p, e = x - p, n = y + p, l = f - p;
    return g > e && (g = r, e = x), n > l && (n = y, l = f), [Math.max(g, Math.min(e, d)), Math.max(n, Math.min(l, s))];
  }
  function z(a, d, s) {
    const u = window.Plotly, [r, x] = T(a, d, s);
    t.met_flux_x[c][a] = r, t.met_flux_y[c][a] = x;
    const y = t.met_ids[a];
    for (let f = 0; f < t.edge_indices.length; f++) {
      const p = t.edge_met_ids[c][f] || [], g = t.edge_flux_x[c][f], e = t.edge_flux_y[c][f];
      for (let n = 0; n < p.length; n++) {
        if (p[n] !== y) continue;
        const l = n * 3 + 1;
        l < g.length && l < e.length && (g[l] = r, e[l] = x);
      }
    }
    u.restyle(i, { x: [t.met_flux_x[c]], y: [t.met_flux_y[c]] }, [t.met_idx]), u.restyle(i, { x: t.edge_flux_x[c], y: t.edge_flux_y[c] }, t.edge_indices);
  }
  function O(a, d) {
    const s = i._fullLayout && i._fullLayout.xaxis, u = i._fullLayout && i._fullLayout.yaxis;
    if (!s || !u) return [0, 0];
    const r = s.range || [0, 1], x = u.range || [0, 1], y = r[1] - r[0], f = x[1] - x[0], p = s._length || 1, g = u._length || 1;
    return [a / p * y, -d / g * f];
  }
  function X() {
    const a = window.Plotly, d = t.rxn_indices.map((s, u) => u === c);
    a.restyle(i, { visible: d }, t.rxn_indices), a.restyle(i, {
      x: t.edge_flux_x[c],
      y: t.edge_flux_y[c]
    }, t.edge_indices), a.restyle(i, {
      x: [t.met_flux_x[c]],
      y: [t.met_flux_y[c]]
    }, [t.met_idx]);
  }
  function Y() {
    const a = window.Plotly;
    i.on("plotly_buttonclicked", (e) => {
      if (!e || !e.button || !e.button.label) return;
      const n = e.button.label;
      if (n === "Enable drag") {
        w = !0, i.style.cursor = "crosshair";
        return;
      }
      if (n === "Disable drag") {
        w = !1, L = null, _ = null, i.style.cursor = "";
        return;
      }
      const l = t.view_labels.indexOf(n);
      l >= 0 && (c = l, X());
    }), i.on("plotly_hover", (e) => {
      if (L = null, !w || !e || !e.points || !e.points.length) return;
      const n = e.points[0];
      n.curveNumber === t.met_idx && Number.isInteger(n.pointNumber) && (L = n.pointNumber, i.style.cursor = "grab");
    }), i.on("plotly_unhover", () => {
      L = null, Number.isInteger(_) || (i.style.cursor = w ? "crosshair" : "");
    }), i.addEventListener("mousedown", (e) => {
      !w || e.button !== 0 || !e.shiftKey || Number.isInteger(L) && (_ = L, v = e.clientX, M = e.clientY, i.style.cursor = "grabbing", e.preventDefault(), e.stopPropagation(), typeof e.stopImmediatePropagation == "function" && e.stopImmediatePropagation());
    }, !0), window.addEventListener("mousemove", (e) => {
      if (!Number.isInteger(_)) return;
      if (v === null || M === null) {
        v = e.clientX, M = e.clientY;
        return;
      }
      const n = e.clientX - v, l = e.clientY - M;
      v = e.clientX, M = e.clientY;
      const [m, b] = O(n, l), o = t.met_flux_x[c][_], h = t.met_flux_y[c][_];
      z(_, o + m, h + b);
    }), window.addEventListener("mouseup", () => {
      Number.isInteger(_) && (_ = null, v = null, M = null, i.style.cursor = w ? "crosshair" : "");
    }), w && (i.style.cursor = "crosshair");
    const d = document.createElement("div");
    d.style.cssText = "position:absolute;top:10px;right:240px;z-index:1000;background:white;border:1px solid #aab7b8;border-radius:4px;padding:4px 8px;box-shadow:0 2px 6px rgba(0,0,0,0.18);width:230px;font-family:Arial,sans-serif;", d.innerHTML = '<input type="text" placeholder="🔍 Search reaction / metabolite…" style="width:100%;border:none;outline:none;font-size:11px;padding:2px 0;box-sizing:border-box;"><div style="display:none;max-height:190px;overflow-y:auto;margin-top:3px;"></div>';
    const s = i.parentElement || document.body;
    s.style.position = "relative", s.appendChild(d);
    const u = d.querySelector("input"), r = d.querySelector("div");
    function x() {
      a.restyle(i, { x: [[]], y: [[]] }, [t.highlight_idx]);
    }
    function y(e, n) {
      a.restyle(i, { x: [e], y: [n] }, [t.highlight_idx]);
    }
    function f(e, n) {
      if (!e.length) return;
      const l = Math.min(...e), m = Math.max(...e), b = Math.min(...n), o = Math.max(...n), h = Math.max((m - l + o - b) * 0.25, 2);
      a.relayout(i, {
        "xaxis.range": [l - h, m + h],
        "yaxis.range": [b - h, o + h]
      });
    }
    function p(e, n) {
      const l = e === "rxn" ? t.rxn_x[n] : t.met_flux_x[c][n], m = e === "rxn" ? t.rxn_y[n] : t.met_flux_y[c][n];
      y([l], [m]), f([l], [m]), r.style.display = "none", u.value = e === "rxn" ? t.rxn_names[n] || t.rxn_ids[n] : t.met_names[n] || t.met_ids[n];
    }
    function g(e) {
      if (r.innerHTML = "", !e.trim()) {
        x(), r.style.display = "none";
        return;
      }
      const n = e.toLowerCase(), l = [];
      for (let o = 0; o < t.rxn_ids.length && l.length < 20; o++)
        (t.rxn_ids[o].toLowerCase().indexOf(n) >= 0 || t.rxn_names[o].toLowerCase().indexOf(n) >= 0) && l.push({ type: "rxn", idx: o });
      for (let o = 0; o < t.met_ids.length && l.length < 20; o++)
        (t.met_ids[o].toLowerCase().indexOf(n) >= 0 || t.met_names[o].toLowerCase().indexOf(n) >= 0) && l.push({ type: "met", idx: o });
      if (!l.length) {
        r.innerHTML = '<div style="padding:4px 8px;color:#999;font-size:11px;">No matches</div>', r.style.display = "block", x();
        return;
      }
      const m = [], b = [];
      for (const o of l)
        m.push(o.type === "rxn" ? t.rxn_x[o.idx] : t.met_flux_x[c][o.idx]), b.push(o.type === "rxn" ? t.rxn_y[o.idx] : t.met_flux_y[c][o.idx]);
      y(m, b), l.length === 1 && f(m, b);
      for (const o of l) {
        const h = o.type === "rxn" ? t.rxn_ids[o.idx] : t.met_ids[o.idx], E = o.type === "rxn" ? t.rxn_names[o.idx] : t.met_names[o.idx], B = o.type === "rxn" ? "rxn" : "met", P = document.createElement("div");
        P.style.cssText = "padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #f0f0f0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;", P.title = h + (E && E !== h ? " — " + E : ""), P.innerHTML = "<b>" + h + "</b>" + (E && E !== h ? ' <span style="color:#555;">— ' + E + "</span>" : "") + ' <span style="color:#aaa;font-size:10px;">[' + B + "]</span>", P.onmouseenter = function() {
          this.style.background = "#f4f6f7";
        }, P.onmouseleave = function() {
          this.style.background = "";
        }, P.onclick = /* @__PURE__ */ ((R, j) => () => p(R, j))(o.type, o.idx), r.appendChild(P);
      }
      r.style.display = "block";
    }
    u.addEventListener("input", function() {
      g(this.value);
    }), u.addEventListener("keydown", function(e) {
      if (e.key === "Escape")
        this.value = "", x(), r.style.display = "none";
      else if (e.key === "Enter") {
        const n = r.querySelector("div[style*=cursor]");
        n && n.click();
      }
    }), document.addEventListener("click", (e) => {
      d.contains(e.target) || (r.style.display = "none");
    });
  }
  function D() {
    const a = window.Plotly, d = S.get("_figure_json");
    if (!d || d === "{}") return;
    let s;
    try {
      s = JSON.parse(d);
    } catch {
      return;
    }
    t = s._interactivity, delete s._interactivity;
    const u = s.data || [], r = s.layout || {}, x = I.getBoundingClientRect().width || 900;
    r.width = x, r.height = Math.max(r.height || 940, 600);
    const y = {
      scrollZoom: !0,
      displaylogo: !1,
      toImageButtonOptions: { format: "svg", filename: "flux_network" }
    };
    i ? a.react(i, u, r, y) : a.newPlot(k, u, r, y).then(() => {
      i = k, t && (c = 0, Y());
    });
  }
  N().then(() => {
    C = !0, D();
  }), S.on("change:_figure_json", () => {
    C && D();
  }), new ResizeObserver((a) => {
    for (const d of a) {
      const s = d.contentRect.width;
      s > 0 && i && C && window.Plotly.relayout(i, { width: s });
    }
  }).observe(I);
}
const A = { render: H };
export {
  A as default
};
