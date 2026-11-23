"use strict";

// The dbSNP lane tab: branch table, per-branch turn charts against the SNPlocs baseline, turn log.
(function () {
  const root = document.getElementById("dbsnp");
  const fmt = new Intl.NumberFormat("en-GB", {maximumFractionDigits: 2});
  const when = new Intl.DateTimeFormat("en-GB", {dateStyle: "medium", timeStyle: "short", timeZone: "UTC"});
  const show = v => (v === null || v === undefined ? "–" : typeof v === "number" ? fmt.format(v) : v);
  const utc = s => {
    if (!s) return "–";
    const m = /^(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})(\d{2})Z$/.exec(s);
    return when.format(m ? new Date(Date.UTC(+m[1], m[2] - 1, +m[3], +m[4], +m[5], +m[6])) : new Date(s));
  };

  function cell(row, text, cls) {
    const td = document.createElement("td");
    td.textContent = show(text);
    if (cls) td.className = cls;
    row.append(td);
    return td;
  }

  function svgEl(parent, name, attrs, text) {
    const e = document.createElementNS("http://www.w3.org/2000/svg", name);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (text !== undefined) e.textContent = text;
    parent.append(e);
    return e;
  }

  // Eight id-lookup branches in fixed palette order (validated on white: dataviz validate_palette.js).
  // Colour follows the branch, never its rank; a marker shape repeats the identity for CVD readers.
  const SERIES = [
    ["canonical", "#2a78d6", "circle"], ["duckdb-parquet", "#eb6834", "square"],
    ["duckdb-native", "#1baf7a", "diamond"], ["arrow-parquet", "#eda100", "triangle"],
    ["flat-mmap", "#e87ba4", "circle"], ["hdf5", "#008300", "square"],
    ["vortex", "#6250d6", "diamond"], ["hail", "#e34948", "triangle"],
  ];

  function marker(parent, shape, x, y, colour, filled) {
    const paint = {fill: filled ? colour : "white", stroke: colour, "stroke-width": 2};
    const r = 4.5;
    if (shape === "square") return svgEl(parent, "rect", {...paint, x: x - r * .85, y: y - r * .85, width: r * 1.7, height: r * 1.7});
    if (shape === "diamond") return svgEl(parent, "polygon", {...paint, points: `${x},${y - r * 1.2} ${x + r * 1.2},${y} ${x},${y + r * 1.2} ${x - r * 1.2},${y}`});
    if (shape === "triangle") return svgEl(parent, "polygon", {...paint, points: `${x},${y - r * 1.2} ${x + r * 1.1},${y + r * .8} ${x - r * 1.1},${y + r * .8}`});
    return svgEl(parent, "circle", {...paint, cx: x, cy: y, r});
  }

  function chart(svg, series, metric, unit, reference, select) {
    svg.replaceChildren();
    const width = svg.clientWidth, height = svg.clientHeight;
    if (!width) return;
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    const left = 53, right = width - 56, top = 20, bottom = height - 34;
    const values = series.flatMap(s => s.turns.map(t => t[metric])).filter(v => v !== null);
    if (reference !== null && reference !== undefined) values.push(reference);
    if (!values.length) {
      svgEl(svg, "text", {x: width / 2, y: height / 2, "text-anchor": "middle"}, "No turns yet");
      return;
    }
    const n = Math.max(1, ...series.map(s => s.turns.length));
    const raw = Math.max(...values);
    const step = 10 ** Math.floor(Math.log10(raw / 4));
    const interval = [1, 2, 2.5, 5, 10].map(v => v * step).find(v => v >= raw / 4);
    const max = interval * Math.ceil(raw * 1.05 / interval);
    const y = v => bottom - (v / max) * (bottom - top);
    const x = i => left + (n === 1 ? 0.5 : i / (n - 1)) * (right - left);
    for (let v = 0; v <= max + interval / 10; v += interval) {
      svgEl(svg, "line", {x1: left, x2: right, y1: y(v), y2: y(v), stroke: "#e1e8e2"});
      svgEl(svg, "text", {x: right + 7, y: y(v) + 4}, fmt.format(v));
    }
    svgEl(svg, "text", {x: right + 7, y: 12}, unit);
    const ticks = Math.min(n, width < 500 ? 4 : 8);
    for (let t = 0; t < ticks; t++) {
      const i = Math.round(t * (n - 1) / Math.max(ticks - 1, 1));
      svgEl(svg, "text", {x: x(i), y: bottom + 16, "text-anchor": "middle"}, i + 1);
    }
    svgEl(svg, "text", {x: (left + right) / 2, y: height - 3, "text-anchor": "middle"}, "Turn within branch");
    if (reference !== null && reference !== undefined) {
      svgEl(svg, "line", {x1: left, x2: right, y1: y(reference), y2: y(reference), stroke: "#405e80", "stroke-width": 2, "stroke-dasharray": "5 4"});
      svgEl(svg, "text", {x: left - 6, y: y(reference) + 4, "text-anchor": "end"}, "SNPlocs");
    }
    for (const s of series) {
      const pts = s.turns.map((t, i) => [x(i), t[metric] === null ? null : y(t[metric])]).filter(p => p[1] !== null);
      if (pts.length > 1) svgEl(svg, "path", {d: "M " + pts.map(p => p.join(" ")).join(" L "), fill: "none", stroke: s.colour, "stroke-width": 2, opacity: .85});
    }
    for (const s of series) {
      s.turns.forEach((t, i) => {
        if (t[metric] === null) return;
        const g = svgEl(svg, "g", {class: "point", tabindex: 0, role: "button",
          "aria-label": `${s.branch} turn ${i + 1}, ${fmt.format(t[metric])} ${unit}, ${t.verdict || ""}`});
        marker(g, s.shape, x(i), y(t[metric]), s.colour, t.verdict === "better" || i === 0);
        svgEl(g, "circle", {cx: x(i), cy: y(t[metric]), r: 10, fill: "transparent"});  // hit target larger than the mark
        const pick = () => select(s, t, i);
        g.addEventListener("pointerenter", pick); g.addEventListener("focus", pick); g.addEventListener("click", pick);
      });
    }
  }

  async function start() {
    const meta = root.querySelector("[data-dbsnp-meta]");
    let data;
    try {
      const response = await fetch("dbsnp.json?v=1", {cache: "no-cache"});
      if (!response.ok) throw new Error("no lane data");
      data = await response.json();
    } catch (error) {
      meta.textContent = "Lane data is not published yet.";
      return;
    }
    const b = data.baseline;
    meta.textContent = `dbSNP ${data.release} · ${data.input} · updated ${utc(data.generated)} UTC` +
      (b ? ` · SNPlocs ${fmt.format(b.cold_s)} s cold, ${fmt.format(b.cold_gb)} GB` : " · SNPlocs baseline not measured yet");

    const body = root.querySelector("[data-dbsnp-branches] tbody");
    for (const br of data.branches) {
      const row = document.createElement("tr");
      cell(row, br.branch);
      cell(row, br.description);
      const stage = cell(row, br.stage);
      if (br.stage === "Waiting" && br.entryFailed.length) stage.title = "Entry gate: " + br.entryFailed.join(", ");
      cell(row, br.turns, "num"); cell(row, br.better, "num");
      cell(row, br.bestColdSeconds, "num"); cell(row, br.bestColdGB, "num");
      cell(row, br.gateA, br.gateA ? "verdict-" + br.gateA : "");
      body.append(row);
    }

    const names = Object.fromEntries(data.branches.map(br => [br.backend, br.branch]));
    const series = SERIES.map(([backend, colour, shape]) => ({
      backend, colour, shape, branch: names[backend] || backend,
      turns: data.turns.filter(t => t.backend === backend && t.kind === "screen"),
    }));
    const legend = root.querySelector("[data-dbsnp-legend]");
    for (const s of series) {
      const item = document.createElement("span");
      item.className = `series ${s.shape}`;
      item.style.setProperty("--c", s.colour);
      item.textContent = `${s.branch} (${s.turns.length})`;
      legend.append(item);
    }
    const detail = root.querySelector("[data-dbsnp-detail]");
    const select = (s, t, i) => {
      detail.textContent = `${s.branch} · turn ${i + 1} · ${show(t.cold_s)} s cold · ${show(t.cold_gb)} GB · ` +
        `${show(t.time_ratio)}× time, ${show(t.ram_ratio)}× RAM vs ${t.ref || "reference"} · ${t.verdict || ""}` +
        (t.note ? ` · ${t.note}` : "") + ` · ${utc(t.utc)} UTC`;
    };
    const runtime = root.querySelector("[data-dbsnp-runtime]"), ram = root.querySelector("[data-dbsnp-ram]");
    function paint() {
      chart(runtime, series, "cold_s", "s", b ? b.cold_s : null, select);
      chart(ram, series, "cold_gb", "GB", b ? b.cold_gb : null, select);
    }
    new ResizeObserver(() => requestAnimationFrame(paint)).observe(runtime);
    paint();

    const log = root.querySelector("[data-dbsnp-turns] tbody");
    for (const t of [...data.turns].reverse().slice(0, 40)) {
      const row = document.createElement("tr");
      cell(row, utc(t.utc)); cell(row, names[t.backend] || t.backend); cell(row, t.kind);
      const out = cell(row, t.b0); out.title = t.b0Detail || "";
      cell(row, t.cold_s, "num"); cell(row, t.cold_gb, "num");
      cell(row, t.time_ratio, "num"); cell(row, t.ram_ratio, "num");
      cell(row, t.verdict, t.verdict ? "verdict-" + t.verdict : "");
      cell(row, t.note, "small");
      log.append(row);
    }
    if (!data.turns.length) {
      const row = document.createElement("tr");
      const td = cell(row, "No turns yet. The first is the SNPlocs baseline screen.");
      td.colSpan = 10;
      log.append(row);
    }
  }

  start();
})();
