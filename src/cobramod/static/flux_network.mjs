function B({ model: N, el: k }) {
  const P = document.createElement("div");
  k.appendChild(P);
  let C = !1, i = null, e = null, u = 0, M = null, b = null, v = null, L = null;
  function z() {
    return window.Plotly ? Promise.resolve() : new Promise((d, a) => {
      const s = document.createElement("script");
      s.src = "https://cdn.plot.ly/plotly-2.35.2.min.js", s.onload = () => d(), s.onerror = () => a(new Error("Failed to load Plotly.js")), document.head.appendChild(s);
    });
  }
  function T(d, a, s) {
    const c = e.met_bounds[d] || [-1e9, 1e9, -1e9, 1e9], [r, x, y, p] = c, m = 0.2;
    let h = r + m, t = x - m, n = y + m, l = p - m;
    return h > t && (h = r, t = x), n > l && (n = y, l = p), [Math.max(h, Math.min(t, a)), Math.max(n, Math.min(l, s))];
  }
  function I(d, a, s) {
    const c = window.Plotly, [r, x] = T(d, a, s);
    e.met_flux_x[u][d] = r, e.met_flux_y[u][d] = x;
    const y = e.met_ids[d];
    for (let p = 0; p < e.edge_indices.length; p++) {
      const m = e.edge_met_ids[u][p] || [], h = e.edge_flux_x[u][p], t = e.edge_flux_y[u][p];
      for (let n = 0; n < m.length; n++) {
        if (m[n] !== y)
          continue;
        const l = n * 3 + 1;
        l < h.length && l < t.length && (h[l] = r, t[l] = x);
      }
    }
    c.restyle(i, { x: [e.met_flux_x[u]], y: [e.met_flux_y[u]] }, [e.met_idx]), c.restyle(i, { x: e.edge_flux_x[u], y: e.edge_flux_y[u] }, e.edge_indices);
  }
  function H(d, a) {
    const s = i._fullLayout && i._fullLayout.xaxis, c = i._fullLayout && i._fullLayout.yaxis;
    if (!s || !c)
      return [0, 0];
    const r = s.range || [0, 1], x = c.range || [0, 1], y = r[1] - r[0], p = x[1] - x[0], m = s._length || 1, h = c._length || 1;
    return [d / m * y, -a / h * p];
  }
  function X() {
    const d = window.Plotly, a = e.rxn_indices.map((s, c) => c === u);
    d.restyle(i, { visible: a }, e.rxn_indices), d.restyle(i, {
      x: e.edge_flux_x[u],
      y: e.edge_flux_y[u]
    }, e.edge_indices), d.restyle(i, {
      x: [e.met_flux_x[u]],
      y: [e.met_flux_y[u]]
    }, [e.met_idx]);
  }
  function Y() {
    const d = window.Plotly;
    i.on("plotly_buttonclicked", (t) => {
      if (!t || !t.button || !t.button.label)
        return;
      const n = t.button.label, l = e.view_labels.indexOf(n);
      l >= 0 && (u = l, X());
    }), i.on("plotly_hover", (t) => {
      if (M = null, !t || !t.points || !t.points.length)
        return;
      const n = t.points[0];
      n.curveNumber === e.met_idx && Number.isInteger(n.pointNumber) && (M = n.pointNumber, i.style.cursor = "grab");
    }), i.on("plotly_unhover", () => {
      M = null, Number.isInteger(b) || (i.style.cursor = "crosshair");
    }), i.addEventListener("mousedown", (t) => {
      t.button !== 0 || !t.shiftKey || Number.isInteger(M) && (b = M, v = t.clientX, L = t.clientY, i.style.cursor = "grabbing", t.preventDefault(), t.stopPropagation(), typeof t.stopImmediatePropagation == "function" && t.stopImmediatePropagation());
    }, !0), window.addEventListener("mousemove", (t) => {
      if (!Number.isInteger(b))
        return;
      if (v === null || L === null) {
        v = t.clientX, L = t.clientY;
        return;
      }
      const n = t.clientX - v, l = t.clientY - L;
      v = t.clientX, L = t.clientY;
      const [f, g] = H(n, l), o = e.met_flux_x[u][b], _ = e.met_flux_y[u][b];
      I(b, o + f, _ + g);
    }), window.addEventListener("mouseup", () => {
      Number.isInteger(b) && (b = null, v = null, L = null, i.style.cursor = "crosshair");
    }), i.style.cursor = "crosshair";
    const a = document.createElement("div");
    a.style.cssText = "position:absolute;top:10px;right:240px;z-index:1000;background:white;border:1px solid #aab7b8;border-radius:4px;padding:4px 8px;box-shadow:0 2px 6px rgba(0,0,0,0.18);width:230px;font-family:Arial,sans-serif;", a.innerHTML = '<input type="text" placeholder="🔍 Search reaction / metabolite…" style="width:100%;border:none;outline:none;font-size:11px;padding:2px 0;box-sizing:border-box;"><div style="display:none;max-height:190px;overflow-y:auto;margin-top:3px;"></div>';
    const s = i.parentElement || document.body;
    s.style.position = "relative", s.appendChild(a);
    const c = a.querySelector("input"), r = a.querySelector("div");
    function x() {
      d.restyle(i, { x: [[]], y: [[]] }, [e.highlight_idx]);
    }
    function y(t, n) {
      d.restyle(i, { x: [t], y: [n] }, [e.highlight_idx]);
    }
    function p(t, n) {
      if (!t.length)
        return;
      const l = Math.min(...t), f = Math.max(...t), g = Math.min(...n), o = Math.max(...n), _ = Math.max((f - l + o - g) * 0.25, 2);
      d.relayout(i, {
        "xaxis.range": [l - _, f + _],
        "yaxis.range": [g - _, o + _]
      });
    }
    function m(t, n) {
      const l = t === "rxn" ? e.rxn_x[n] : e.met_flux_x[u][n], f = t === "rxn" ? e.rxn_y[n] : e.met_flux_y[u][n];
      y([l], [f]), p([l], [f]), r.style.display = "none", c.value = t === "rxn" ? e.rxn_names[n] || e.rxn_ids[n] : e.met_names[n] || e.met_ids[n];
    }
    function h(t) {
      if (r.innerHTML = "", !t.trim()) {
        x(), r.style.display = "none";
        return;
      }
      const n = t.toLowerCase(), l = [], f = [], g = [];
      for (let o = 0; o < e.rxn_ids.length; o++)
        (e.rxn_ids[o].toLowerCase().indexOf(n) >= 0 || e.rxn_names[o].toLowerCase().indexOf(n) >= 0) && (l.push(e.rxn_x[o]), f.push(e.rxn_y[o]), g.length < 20 && g.push({ type: "rxn", idx: o }));
      for (let o = 0; o < e.met_ids.length; o++)
        (e.met_ids[o].toLowerCase().indexOf(n) >= 0 || e.met_names[o].toLowerCase().indexOf(n) >= 0) && (l.push(e.met_flux_x[u][o]), f.push(e.met_flux_y[u][o]), g.length < 20 && g.push({ type: "met", idx: o }));
      if (!l.length) {
        r.innerHTML = '<div style="padding:4px 8px;color:#999;font-size:11px;">No matches</div>', r.style.display = "block", x();
        return;
      }
      y(l, f), l.length === 1 && p(l, f);
      for (const o of g) {
        const _ = o.type === "rxn" ? e.rxn_ids[o.idx] : e.met_ids[o.idx], E = o.type === "rxn" ? e.rxn_names[o.idx] : e.met_names[o.idx], j = o.type === "rxn" ? "rxn" : "met", w = document.createElement("div");
        w.style.cssText = "padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #f0f0f0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;", w.title = _ + (E && E !== _ ? " — " + E : ""), w.innerHTML = "<b>" + _ + "</b>" + (E && E !== _ ? ' <span style="color:#555;">— ' + E + "</span>" : "") + ' <span style="color:#aaa;font-size:10px;">[' + j + "]</span>", w.onmouseenter = function() {
          this.style.background = "#f4f6f7";
        }, w.onmouseleave = function() {
          this.style.background = "";
        }, w.onclick = /* @__PURE__ */ ((S, q) => () => m(S, q))(o.type, o.idx), r.appendChild(w);
      }
      r.style.display = "block";
    }
    c.addEventListener("input", function() {
      h(this.value);
    }), c.addEventListener("keydown", function(t) {
      if (t.key === "Escape")
        this.value = "", x(), r.style.display = "none";
      else if (t.key === "Enter") {
        const n = r.querySelector("div[style*=cursor]");
        n && n.click();
      }
    }), document.addEventListener("click", (t) => {
      a.contains(t.target) || (r.style.display = "none");
    });
  }
  function O() {
    const d = window.Plotly, a = N.get("_figure_json");
    if (!a || a === "{}")
      return;
    let s;
    try {
      s = JSON.parse(a);
    } catch {
      return;
    }
    e = s._interactivity, delete s._interactivity;
    const c = s.data || [], r = s.layout || {}, x = k.getBoundingClientRect().width || 900;
    r.width = x, r.height = Math.max(r.height || 940, 600);
    const y = {
      scrollZoom: !0,
      displaylogo: !1,
      toImageButtonOptions: { format: "svg", filename: "flux_network" }
    };
    i ? d.react(i, c, r, y) : d.newPlot(P, c, r, y).then(() => {
      i = P, e && (u = 0, Y());
    });
  }
  z().then(() => {
    C = !0, O();
  }), N.on("change:_figure_json", () => {
    C && O();
  }), new ResizeObserver((d) => {
    for (const a of d) {
      const s = a.contentRect.width;
      s > 0 && i && C && window.Plotly.relayout(i, { width: s });
    }
  }).observe(k);
}
const R = { render: B };
export {
  R as default
};
