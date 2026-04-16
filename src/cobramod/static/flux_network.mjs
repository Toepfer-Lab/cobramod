function D({ model: O, el: P }) {
  const C = document.createElement("div");
  P.appendChild(C);
  let N = !1, i = null, e = null, d = 0, w = !0, L = null, b = null, E = null, M = null;
  function T() {
    return window.Plotly ? Promise.resolve() : new Promise((u, a) => {
      const r = document.createElement("script");
      r.src = "https://cdn.plot.ly/plotly-2.35.2.min.js", r.onload = () => u(), r.onerror = () => a(new Error("Failed to load Plotly.js")), document.head.appendChild(r);
    });
  }
  function I(u, a, r) {
    const c = e.met_bounds[u] || [-1e9, 1e9, -1e9, 1e9], [s, x, p, y] = c, m = 0.2;
    let h = s + m, t = x - m, n = p + m, l = y - m;
    return h > t && (h = s, t = x), n > l && (n = p, l = y), [Math.max(h, Math.min(t, a)), Math.max(n, Math.min(l, r))];
  }
  function H(u, a, r) {
    const c = window.Plotly, [s, x] = I(u, a, r);
    e.met_flux_x[d][u] = s, e.met_flux_y[d][u] = x;
    const p = e.met_ids[u];
    for (let y = 0; y < e.edge_indices.length; y++) {
      const m = e.edge_met_ids[d][y] || [], h = e.edge_flux_x[d][y], t = e.edge_flux_y[d][y];
      for (let n = 0; n < m.length; n++) {
        if (m[n] !== p)
          continue;
        const l = n * 3 + 1;
        l < h.length && l < t.length && (h[l] = s, t[l] = x);
      }
    }
    c.restyle(i, { x: [e.met_flux_x[d]], y: [e.met_flux_y[d]] }, [e.met_idx]), c.restyle(i, { x: e.edge_flux_x[d], y: e.edge_flux_y[d] }, e.edge_indices);
  }
  function X(u, a) {
    const r = i._fullLayout && i._fullLayout.xaxis, c = i._fullLayout && i._fullLayout.yaxis;
    if (!r || !c)
      return [0, 0];
    const s = r.range || [0, 1], x = c.range || [0, 1], p = s[1] - s[0], y = x[1] - x[0], m = r._length || 1, h = c._length || 1;
    return [u / m * p, -a / h * y];
  }
  function Y() {
    const u = window.Plotly, a = e.rxn_indices.map((r, c) => c === d);
    u.restyle(i, { visible: a }, e.rxn_indices), u.restyle(i, {
      x: e.edge_flux_x[d],
      y: e.edge_flux_y[d]
    }, e.edge_indices), u.restyle(i, {
      x: [e.met_flux_x[d]],
      y: [e.met_flux_y[d]]
    }, [e.met_idx]);
  }
  function j() {
    const u = window.Plotly;
    i.on("plotly_buttonclicked", (t) => {
      if (!t || !t.button || !t.button.label)
        return;
      const n = t.button.label;
      if (n === "Enable drag") {
        w = !0, i.style.cursor = "crosshair";
        return;
      }
      if (n === "Disable drag") {
        w = !1, L = null, b = null, i.style.cursor = "";
        return;
      }
      const l = e.view_labels.indexOf(n);
      l >= 0 && (d = l, Y());
    }), i.on("plotly_hover", (t) => {
      if (L = null, !w || !t || !t.points || !t.points.length)
        return;
      const n = t.points[0];
      n.curveNumber === e.met_idx && Number.isInteger(n.pointNumber) && (L = n.pointNumber, i.style.cursor = "grab");
    }), i.on("plotly_unhover", () => {
      L = null, Number.isInteger(b) || (i.style.cursor = w ? "crosshair" : "");
    }), i.addEventListener("mousedown", (t) => {
      !w || t.button !== 0 || !t.shiftKey || Number.isInteger(L) && (b = L, E = t.clientX, M = t.clientY, i.style.cursor = "grabbing", t.preventDefault(), t.stopPropagation(), typeof t.stopImmediatePropagation == "function" && t.stopImmediatePropagation());
    }, !0), window.addEventListener("mousemove", (t) => {
      if (!Number.isInteger(b))
        return;
      if (E === null || M === null) {
        E = t.clientX, M = t.clientY;
        return;
      }
      const n = t.clientX - E, l = t.clientY - M;
      E = t.clientX, M = t.clientY;
      const [f, g] = X(n, l), o = e.met_flux_x[d][b], _ = e.met_flux_y[d][b];
      H(b, o + f, _ + g);
    }), window.addEventListener("mouseup", () => {
      Number.isInteger(b) && (b = null, E = null, M = null, i.style.cursor = w ? "crosshair" : "");
    }), w && (i.style.cursor = "crosshair");
    const a = document.createElement("div");
    a.style.cssText = "position:absolute;top:10px;right:240px;z-index:1000;background:white;border:1px solid #aab7b8;border-radius:4px;padding:4px 8px;box-shadow:0 2px 6px rgba(0,0,0,0.18);width:230px;font-family:Arial,sans-serif;", a.innerHTML = '<input type="text" placeholder="🔍 Search reaction / metabolite…" style="width:100%;border:none;outline:none;font-size:11px;padding:2px 0;box-sizing:border-box;"><div style="display:none;max-height:190px;overflow-y:auto;margin-top:3px;"></div>';
    const r = i.parentElement || document.body;
    r.style.position = "relative", r.appendChild(a);
    const c = a.querySelector("input"), s = a.querySelector("div");
    function x() {
      u.restyle(i, { x: [[]], y: [[]] }, [e.highlight_idx]);
    }
    function p(t, n) {
      u.restyle(i, { x: [t], y: [n] }, [e.highlight_idx]);
    }
    function y(t, n) {
      if (!t.length)
        return;
      const l = Math.min(...t), f = Math.max(...t), g = Math.min(...n), o = Math.max(...n), _ = Math.max((f - l + o - g) * 0.25, 2);
      u.relayout(i, {
        "xaxis.range": [l - _, f + _],
        "yaxis.range": [g - _, o + _]
      });
    }
    function m(t, n) {
      const l = t === "rxn" ? e.rxn_x[n] : e.met_flux_x[d][n], f = t === "rxn" ? e.rxn_y[n] : e.met_flux_y[d][n];
      p([l], [f]), y([l], [f]), s.style.display = "none", c.value = t === "rxn" ? e.rxn_names[n] || e.rxn_ids[n] : e.met_names[n] || e.met_ids[n];
    }
    function h(t) {
      if (s.innerHTML = "", !t.trim()) {
        x(), s.style.display = "none";
        return;
      }
      const n = t.toLowerCase(), l = [], f = [], g = [];
      for (let o = 0; o < e.rxn_ids.length; o++)
        (e.rxn_ids[o].toLowerCase().indexOf(n) >= 0 || e.rxn_names[o].toLowerCase().indexOf(n) >= 0) && (l.push(e.rxn_x[o]), f.push(e.rxn_y[o]), g.length < 20 && g.push({ type: "rxn", idx: o }));
      for (let o = 0; o < e.met_ids.length; o++)
        (e.met_ids[o].toLowerCase().indexOf(n) >= 0 || e.met_names[o].toLowerCase().indexOf(n) >= 0) && (l.push(e.met_flux_x[d][o]), f.push(e.met_flux_y[d][o]), g.length < 20 && g.push({ type: "met", idx: o }));
      if (!l.length) {
        s.innerHTML = '<div style="padding:4px 8px;color:#999;font-size:11px;">No matches</div>', s.style.display = "block", x();
        return;
      }
      p(l, f), l.length === 1 && y(l, f);
      for (const o of g) {
        const _ = o.type === "rxn" ? e.rxn_ids[o.idx] : e.met_ids[o.idx], k = o.type === "rxn" ? e.rxn_names[o.idx] : e.met_names[o.idx], S = o.type === "rxn" ? "rxn" : "met", v = document.createElement("div");
        v.style.cssText = "padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #f0f0f0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;", v.title = _ + (k && k !== _ ? " — " + k : ""), v.innerHTML = "<b>" + _ + "</b>" + (k && k !== _ ? ' <span style="color:#555;">— ' + k + "</span>" : "") + ' <span style="color:#aaa;font-size:10px;">[' + S + "]</span>", v.onmouseenter = function() {
          this.style.background = "#f4f6f7";
        }, v.onmouseleave = function() {
          this.style.background = "";
        }, v.onclick = /* @__PURE__ */ ((q, B) => () => m(q, B))(o.type, o.idx), s.appendChild(v);
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
  function z() {
    const u = window.Plotly, a = O.get("_figure_json");
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
    const p = {
      scrollZoom: !0,
      displaylogo: !1,
      toImageButtonOptions: { format: "svg", filename: "flux_network" }
    };
    i ? u.react(i, c, s, p) : u.newPlot(C, c, s, p).then(() => {
      i = C, e && (d = 0, j());
    });
  }
  T().then(() => {
    N = !0, z();
  }), O.on("change:_figure_json", () => {
    N && z();
  }), new ResizeObserver((u) => {
    for (const a of u) {
      const r = a.contentRect.width;
      r > 0 && i && N && window.Plotly.relayout(i, { width: r });
    }
  }).observe(P);
}
const R = { render: D };
export {
  R as default
};
