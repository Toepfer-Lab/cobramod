function D({ model: I, el: P }) {
  const C = document.createElement("div");
  P.appendChild(C);
  let N = !1, i = null, e = null, u = 0, w = !0, L = null, g = null, E = null, M = null;
  function z() {
    return window.Plotly ? Promise.resolve() : new Promise((d, a) => {
      const r = document.createElement("script");
      r.src = "https://cdn.plot.ly/plotly-2.35.2.min.js", r.onload = () => d(), r.onerror = () => a(new Error("Failed to load Plotly.js")), document.head.appendChild(r);
    });
  }
  function T(d, a, r) {
    const c = e.met_bounds[d] || [-1e9, 1e9, -1e9, 1e9], [s, x, f, y] = c, p = 0.2;
    let h = s + p, t = x - p, n = f + p, l = y - p;
    return h > t && (h = s, t = x), n > l && (n = f, l = y), [Math.max(h, Math.min(t, a)), Math.max(n, Math.min(l, r))];
  }
  function S(d, a, r) {
    const c = window.Plotly, [s, x] = T(d, a, r);
    e.met_flux_x[u][d] = s, e.met_flux_y[u][d] = x;
    const f = e.met_ids[d];
    for (let y = 0; y < e.edge_indices.length; y++) {
      const p = e.edge_met_ids[u][y] || [], h = e.edge_flux_x[u][y], t = e.edge_flux_y[u][y];
      for (let n = 0; n < p.length; n++) {
        if (p[n] !== f)
          continue;
        const l = n * 3 + 1;
        l < h.length && l < t.length && (h[l] = s, t[l] = x);
      }
    }
    c.restyle(i, { x: [e.met_flux_x[u]], y: [e.met_flux_y[u]] }, [e.met_idx]), c.restyle(i, { x: e.edge_flux_x[u], y: e.edge_flux_y[u] }, e.edge_indices);
  }
  function B(d, a) {
    const r = i._fullLayout && i._fullLayout.xaxis, c = i._fullLayout && i._fullLayout.yaxis;
    if (!r || !c)
      return [0, 0];
    const s = r.range || [0, 1], x = c.range || [0, 1], f = s[1] - s[0], y = x[1] - x[0], p = r._length || 1, h = c._length || 1;
    return [d / p * f, -a / h * y];
  }
  function H() {
    const d = window.Plotly, a = e.rxn_indices.map((r, c) => c === u);
    d.restyle(i, { visible: a }, e.rxn_indices), d.restyle(i, {
      x: e.edge_flux_x[u],
      y: e.edge_flux_y[u]
    }, e.edge_indices), d.restyle(i, {
      x: [e.met_flux_x[u]],
      y: [e.met_flux_y[u]]
    }, [e.met_idx]);
  }
  function R() {
    const d = window.Plotly;
    i.on("plotly_buttonclicked", (t) => {
      if (!t || !t.button || !t.button.label)
        return;
      const n = t.button.label;
      if (n === "Enable drag") {
        w = !0, i.style.cursor = "crosshair";
        return;
      }
      if (n === "Disable drag") {
        w = !1, L = null, g = null, i.style.cursor = "";
        return;
      }
      const l = e.view_labels.indexOf(n);
      l >= 0 && (u = l, H());
    }), i.on("plotly_hover", (t) => {
      if (L = null, !w || !t || !t.points || !t.points.length)
        return;
      const n = t.points[0];
      n.curveNumber === e.met_idx && Number.isInteger(n.pointNumber) && (L = n.pointNumber, i.style.cursor = "grab");
    }), i.on("plotly_unhover", () => {
      L = null, Number.isInteger(g) || (i.style.cursor = w ? "crosshair" : "");
    }), i.addEventListener("mousedown", (t) => {
      !w || t.button !== 0 || !t.shiftKey || Number.isInteger(L) && (g = L, E = t.clientX, M = t.clientY, i.style.cursor = "grabbing", t.preventDefault(), t.stopPropagation(), typeof t.stopImmediatePropagation == "function" && t.stopImmediatePropagation());
    }, !0), window.addEventListener("mousemove", (t) => {
      if (!Number.isInteger(g))
        return;
      if (E === null || M === null) {
        E = t.clientX, M = t.clientY;
        return;
      }
      const n = t.clientX - E, l = t.clientY - M;
      E = t.clientX, M = t.clientY;
      const [m, b] = B(n, l), o = e.met_flux_x[u][g], _ = e.met_flux_y[u][g];
      S(g, o + m, _ + b);
    }), window.addEventListener("mouseup", () => {
      Number.isInteger(g) && (g = null, E = null, M = null, i.style.cursor = w ? "crosshair" : "");
    }), w && (i.style.cursor = "crosshair");
    const a = document.createElement("div");
    a.style.cssText = "position:absolute;top:10px;right:240px;z-index:1000;background:white;border:1px solid #aab7b8;border-radius:4px;padding:4px 8px;box-shadow:0 2px 6px rgba(0,0,0,0.18);width:230px;font-family:Arial,sans-serif;", a.innerHTML = '<input id="fluxSrchIn" type="text" placeholder="🔍 Search reaction / metabolite…" style="width:100%;border:none;outline:none;font-size:11px;padding:2px 0;box-sizing:border-box;"><div id="fluxSrchRes" style="display:none;max-height:190px;overflow-y:auto;margin-top:3px;"></div>';
    const r = i.parentElement || document.body;
    r.style.position = "relative", r.appendChild(a);
    const c = document.getElementById("fluxSrchIn"), s = document.getElementById("fluxSrchRes");
    function x() {
      d.restyle(i, { x: [[]], y: [[]] }, [e.highlight_idx]);
    }
    function f(t, n) {
      d.restyle(i, { x: [t], y: [n] }, [e.highlight_idx]);
    }
    function y(t, n) {
      if (!t.length)
        return;
      const l = Math.min(...t), m = Math.max(...t), b = Math.min(...n), o = Math.max(...n), _ = Math.max((m - l + o - b) * 0.25, 2);
      d.relayout(i, {
        "xaxis.range": [l - _, m + _],
        "yaxis.range": [b - _, o + _]
      });
    }
    function p(t, n) {
      const l = t === "rxn" ? e.rxn_x[n] : e.met_flux_x[u][n], m = t === "rxn" ? e.rxn_y[n] : e.met_flux_y[u][n];
      f([l], [m]), y([l], [m]), s.style.display = "none", c.value = t === "rxn" ? e.rxn_names[n] || e.rxn_ids[n] : e.met_names[n] || e.met_ids[n];
    }
    function h(t) {
      if (s.innerHTML = "", !t.trim()) {
        x(), s.style.display = "none";
        return;
      }
      const n = t.toLowerCase(), l = [];
      for (let o = 0; o < e.rxn_ids.length && l.length < 20; o++)
        (e.rxn_ids[o].toLowerCase().indexOf(n) >= 0 || e.rxn_names[o].toLowerCase().indexOf(n) >= 0) && l.push({ type: "rxn", idx: o });
      for (let o = 0; o < e.met_ids.length && l.length < 20; o++)
        (e.met_ids[o].toLowerCase().indexOf(n) >= 0 || e.met_names[o].toLowerCase().indexOf(n) >= 0) && l.push({ type: "met", idx: o });
      if (!l.length) {
        s.innerHTML = '<div style="padding:4px 8px;color:#999;font-size:11px;">No matches</div>', s.style.display = "block", x();
        return;
      }
      const m = [], b = [];
      for (const o of l)
        m.push(o.type === "rxn" ? e.rxn_x[o.idx] : e.met_flux_x[u][o.idx]), b.push(o.type === "rxn" ? e.rxn_y[o.idx] : e.met_flux_y[u][o.idx]);
      f(m, b), l.length === 1 && y(m, b);
      for (const o of l) {
        const _ = o.type === "rxn" ? e.rxn_ids[o.idx] : e.met_ids[o.idx], k = o.type === "rxn" ? e.rxn_names[o.idx] : e.met_names[o.idx], X = o.type === "rxn" ? "rxn" : "met", v = document.createElement("div");
        v.style.cssText = "padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #f0f0f0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;", v.title = _ + (k && k !== _ ? " — " + k : ""), v.innerHTML = "<b>" + _ + "</b>" + (k && k !== _ ? ' <span style="color:#555;">— ' + k + "</span>" : "") + ' <span style="color:#aaa;font-size:10px;">[' + X + "]</span>", v.onmouseenter = function() {
          this.style.background = "#f4f6f7";
        }, v.onmouseleave = function() {
          this.style.background = "";
        }, v.onclick = /* @__PURE__ */ ((Y, j) => () => p(Y, j))(o.type, o.idx), s.appendChild(v);
      }
      s.style.display = "block";
    }
    c.addEventListener("input", function() {
      h(this.value);
    }), c.addEventListener("keydown", function(t) {
      if (t.key === "Escape")
        this.value = "", x(), s.style.display = "none";
      else if (t.key === "Enter") {
        const n = s.querySelector("div[style*=cursor]");
        n && n.click();
      }
    }), document.addEventListener("click", (t) => {
      a.contains(t.target) || (s.style.display = "none");
    });
  }
  function O() {
    const d = window.Plotly, a = I.get("_figure_json");
    if (!a || a === "{}")
      return;
    let r;
    try {
      r = JSON.parse(a);
    } catch {
      return;
    }
    e = r._interactivity, delete r._interactivity;
    const c = r.data || [], s = r.layout || {}, x = P.getBoundingClientRect().width || 900;
    s.width = x, s.height = Math.max(s.height || 940, 600);
    const f = {
      scrollZoom: !0,
      displaylogo: !1,
      toImageButtonOptions: { format: "svg", filename: "flux_network" }
    };
    i ? d.react(i, c, s, f) : d.newPlot(C, c, s, f).then(() => {
      i = C, e && (u = 0, R());
    });
  }
  z().then(() => {
    N = !0, O();
  }), I.on("change:_figure_json", () => {
    N && O();
  }), new ResizeObserver((d) => {
    for (const a of d) {
      const r = a.contentRect.width;
      r > 0 && i && N && window.Plotly.relayout(i, { width: r });
    }
  }).observe(P);
}
const A = { render: D };
export {
  A as default
};
