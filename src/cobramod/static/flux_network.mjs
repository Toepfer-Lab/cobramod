function S({ model: T, el: P }) {
  const N = document.createElement("div");
  P.appendChild(N);
  let C = !1, i = null, e = null, a = 0, E = null, f = null, M = null, k = null, b = null, w = 0, v = 0;
  function I() {
    return window.Plotly ? Promise.resolve() : new Promise((l, d) => {
      const r = document.createElement("script");
      r.src = "https://cdn.plot.ly/plotly-2.35.2.min.js", r.onload = () => l(), r.onerror = () => d(new Error("Failed to load Plotly.js")), document.head.appendChild(r);
    });
  }
  function H(l, d, r) {
    const m = e.met_bounds[l] || [-1e9, 1e9, -1e9, 1e9], [s, p, y, c] = m, t = 0.2;
    let n = s + t, u = p - t, x = y + t, _ = c - t;
    return n > u && (n = s, u = p), x > _ && (x = y, _ = c), [Math.max(n, Math.min(u, d)), Math.max(x, Math.min(_, r))];
  }
  function j(l, d, r) {
    const m = window.Plotly, [s, p] = H(l, d, r);
    e.met_flux_x[a][l] = s, e.met_flux_y[a][l] = p;
    const y = e.met_ids[l];
    for (let c = 0; c < e.edge_indices.length; c++) {
      const t = e.edge_met_ids[a][c] || [], n = e.edge_flux_x[a][c], u = e.edge_flux_y[a][c];
      for (let x = 0; x < t.length; x++) {
        if (t[x] !== y)
          continue;
        const _ = x * 3 + 1;
        _ < n.length && _ < u.length && (n[_] = s, u[_] = p);
      }
    }
    m.restyle(i, { x: [e.met_flux_x[a], ...e.edge_flux_x[a]], y: [e.met_flux_y[a], ...e.edge_flux_y[a]] }, [e.met_idx, ...e.edge_indices]);
  }
  function q(l, d) {
    const r = i._fullLayout && i._fullLayout.xaxis, m = i._fullLayout && i._fullLayout.yaxis;
    if (!r || !m)
      return [0, 0];
    const s = r.range || [0, 1], p = m.range || [0, 1], y = s[1] - s[0], c = p[1] - p[0], t = r._length || 1, n = m._length || 1;
    return [l / t * y, -d / n * c];
  }
  function A() {
    const l = window.Plotly, d = e.rxn_indices.map((r, m) => m === a);
    l.restyle(i, { visible: d }, e.rxn_indices), l.restyle(i, {
      x: e.edge_flux_x[a],
      y: e.edge_flux_y[a]
    }, e.edge_indices), l.restyle(i, {
      x: [e.met_flux_x[a]],
      y: [e.met_flux_y[a]]
    }, [e.met_idx]);
  }
  function D() {
    const l = window.Plotly;
    i.on("plotly_buttonclicked", (t) => {
      if (!t || !t.button || !t.button.label)
        return;
      const n = e.view_labels.indexOf(t.button.label);
      n >= 0 && (a = n, A());
    }), i.on("plotly_hover", (t) => {
      if (E = null, !t || !t.points || !t.points.length)
        return;
      const n = t.points[0];
      n.curveNumber === e.met_idx && Number.isInteger(n.pointNumber) && (E = n.pointNumber, i.style.cursor = "grab");
    }), i.on("plotly_unhover", () => {
      E = null, Number.isInteger(f) || (i.style.cursor = "crosshair");
    }), i.addEventListener("mousedown", (t) => {
      t.button !== 0 || !t.shiftKey || !Number.isInteger(E) || (f = E, M = t.clientX, k = t.clientY, w = v = 0, i.style.cursor = "grabbing", t.preventDefault(), t.stopPropagation(), typeof t.stopImmediatePropagation == "function" && t.stopImmediatePropagation());
    }, !0), window.addEventListener("mousemove", (t) => {
      Number.isInteger(f) && (w += t.clientX - M, v += t.clientY - k, M = t.clientX, k = t.clientY, b === null && (b = requestAnimationFrame(() => {
        if (b = null, !Number.isInteger(f)) {
          w = v = 0;
          return;
        }
        const [n, u] = q(w, v);
        w = v = 0, j(f, e.met_flux_x[a][f] + n, e.met_flux_y[a][f] + u);
      })));
    }), window.addEventListener("mouseup", () => {
      Number.isInteger(f) && (b !== null && (cancelAnimationFrame(b), b = null), w = v = 0, f = null, M = null, k = null, i.style.cursor = "crosshair");
    }), i.style.cursor = "crosshair";
    const d = document.createElement("div");
    d.style.cssText = "position:absolute;top:10px;right:240px;z-index:1000;background:white;border:1px solid #aab7b8;border-radius:4px;padding:4px 8px;box-shadow:0 2px 6px rgba(0,0,0,0.18);width:230px;font-family:Arial,sans-serif;", d.innerHTML = '<input type="text" placeholder="🔍 Search reaction / metabolite…" style="width:100%;border:none;outline:none;font-size:11px;padding:2px 0;box-sizing:border-box;"><div style="display:none;max-height:190px;overflow-y:auto;margin-top:3px;"></div>';
    const r = i.parentElement || document.body;
    r.style.position = "relative", r.appendChild(d);
    const m = d.querySelector("input"), s = d.querySelector("div");
    function p(t, n) {
      if (!t.length)
        return;
      const u = Math.min(...t), x = Math.max(...t), _ = Math.min(...n), o = Math.max(...n), h = Math.max((x - u + o - _) * 0.25, 2);
      l.relayout(i, {
        "xaxis.range": [u - h, x + h],
        "yaxis.range": [_ - h, o + h]
      });
    }
    function y(t, n) {
      const u = t === "rxn" ? e.rxn_x[n] : e.met_flux_x[a][n], x = t === "rxn" ? e.rxn_y[n] : e.met_flux_y[a][n];
      l.restyle(i, { x: [[u]], y: [[x]] }, [e.highlight_idx]), p([u], [x]), s.style.display = "none", m.value = t === "rxn" ? e.rxn_names[n] || e.rxn_ids[n] : e.met_names[n] || e.met_ids[n];
    }
    function c(t) {
      if (s.innerHTML = "", !t.trim()) {
        l.restyle(i, { x: [[]], y: [[]] }, [e.highlight_idx]), s.style.display = "none";
        return;
      }
      const n = t.toLowerCase(), u = [], x = [], _ = [];
      for (let o = 0; o < e.rxn_ids.length; o++)
        (e._rxn_ids_lc[o].indexOf(n) >= 0 || e._rxn_names_lc[o].indexOf(n) >= 0) && (u.push(e.rxn_x[o]), x.push(e.rxn_y[o]), _.length < 20 && _.push({ type: "rxn", idx: o }));
      for (let o = 0; o < e.met_ids.length; o++)
        (e._met_ids_lc[o].indexOf(n) >= 0 || e._met_names_lc[o].indexOf(n) >= 0) && (u.push(e.met_flux_x[a][o]), x.push(e.met_flux_y[a][o]), _.length < 20 && _.push({ type: "met", idx: o }));
      if (!u.length) {
        s.innerHTML = '<div style="padding:4px 8px;color:#999;font-size:11px;">No matches</div>', s.style.display = "block", l.restyle(i, { x: [[]], y: [[]] }, [e.highlight_idx]);
        return;
      }
      l.restyle(i, { x: [u], y: [x] }, [e.highlight_idx]), u.length === 1 && p(u, x);
      for (const o of _) {
        const h = o.type === "rxn" ? e.rxn_ids[o.idx] : e.met_ids[o.idx], L = o.type === "rxn" ? e.rxn_names[o.idx] : e.met_names[o.idx], g = document.createElement("div");
        g.style.cssText = "padding:3px 8px;cursor:pointer;font-size:11px;border-top:1px solid #f0f0f0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;", g.title = h + (L && L !== h ? " — " + L : ""), g.innerHTML = "<b>" + h + "</b>" + (L && L !== h ? ' <span style="color:#555;">— ' + L + "</span>" : "") + ' <span style="color:#aaa;font-size:10px;">[' + o.type + "]</span>", g.onmouseenter = function() {
          this.style.background = "#f4f6f7";
        }, g.onmouseleave = function() {
          this.style.background = "";
        }, g.onclick = () => y(o.type, o.idx), s.appendChild(g);
      }
      s.style.display = "block";
    }
    m.addEventListener("input", function() {
      c(this.value);
    }), m.addEventListener("keydown", function(t) {
      if (t.key === "Escape")
        this.value = "", l.restyle(i, { x: [[]], y: [[]] }, [e.highlight_idx]), s.style.display = "none";
      else if (t.key === "Enter") {
        const n = s.querySelector("div[style*=cursor]");
        n && n.click();
      }
    }), document.addEventListener("click", (t) => {
      d.contains(t.target) || (s.style.display = "none");
    });
  }
  function z() {
    const l = window.Plotly, d = T.get("_figure_json");
    if (!d || d === "{}")
      return;
    let r;
    try {
      r = JSON.parse(d);
    } catch {
      return;
    }
    e = r._interactivity, delete r._interactivity, e._rxn_ids_lc = e.rxn_ids.map((c) => c.toLowerCase()), e._rxn_names_lc = e.rxn_names.map((c) => c.toLowerCase()), e._met_ids_lc = e.met_ids.map((c) => c.toLowerCase()), e._met_names_lc = e.met_names.map((c) => c.toLowerCase());
    const m = r.data || [], s = r.layout || {}, p = P.getBoundingClientRect().width || 900;
    s.width = p, s.height = Math.max(s.height || 940, 600);
    const y = {
      scrollZoom: !0,
      displaylogo: !1,
      toImageButtonOptions: { format: "svg", filename: "flux_network" }
    };
    i ? l.react(i, m, s, y) : l.newPlot(N, m, s, y).then(() => {
      i = N, e && (a = 0, D());
    });
  }
  let O = null;
  I().then(() => {
    C = !0, z();
  }), T.on("change:_figure_json", () => {
    C && z();
  }), new ResizeObserver((l) => {
    for (const d of l) {
      const r = d.contentRect.width;
      r <= 0 || !i || !C || (clearTimeout(O), O = setTimeout(() => window.Plotly.relayout(i, { width: r }), 100));
    }
  }).observe(P);
}
const X = { render: S };
export {
  X as default
};
