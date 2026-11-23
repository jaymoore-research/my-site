"use strict";

// The dbSNP lane tab: branch table, per-branch turn charts against the SNPlocs baseline, turn log.
(function () {
  const root = document.getElementById("dbsnp");
  const fmt = new Intl.NumberFormat("en-GB", {maximumFractionDigits: 2});
  const when = new Intl.DateTimeFormat("en-GB", {dateStyle: "medium", timeStyle: "short", timeZone: "UTC"});
  const verdictColour = {better: "#23744d", same: "#65746e", worse: "#9a3b2e", wrong: "#9a3b2e", baseline: "#405e80"};
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

  function chart(svg, turns, metric, unit, reference) {
    svg.replaceChildren();
    const width = svg.clientWidth, height = svg.clientHeight;
    if (!width) return;
    svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    const left = 53, right = width - 56, top = 20, bottom = height - 34;
    const values = turns.map(t => t[metric]).filter(v => v !== null);
    if (reference !== null && reference !== undefined) values.push(reference);
    if (!values.length) {
      svgEl(svg, "text", {x: width / 2, y: height / 2, "text-anchor": "middle"}, "No turns yet");
      return;
    }
    const max = Math.max(...values) * 1.1;
    const y = v => bottom - (v / max) * (bottom - top);
    const x = i => left + (turns.length === 1 ? 0.5 : i / (turns.length - 1)) * (right - left);
    for (let k = 0; k <= 4; k++) {
      const v = (max * k) / 4;
      svgEl(svg, "line", {x1: left, x2: right, y1: y(v), y2: y(v), stroke: "#e1e8e2"});
      svgEl(svg, "text", {x: right + 7, y: y(v) + 4}, fmt.format(v));
    }
    svgEl(svg, "text", {x: right + 7, y: 12}, unit);
    if (reference !== null && reference !== undefined) {
      svgEl(svg, "line", {x1: left, x2: right, y1: y(reference), y2: y(reference), stroke: "#405e80", "stroke-dasharray": "4 4"});
      svgEl(svg, "text", {x: left - 6, y: y(reference) + 4, "text-anchor": "end"}, "SNPlocs");
    }
    turns.forEach((t, i) => {
      svgEl(svg, "text", {x: x(i), y: bottom + 16, "text-anchor": "middle"}, i + 1);
      if (t[metric] === null) return;
      const dot = svgEl(svg, "circle", {cx: x(i), cy: y(t[metric]), r: 4, fill: verdictColour[t.verdict] || "#65746e"});
      svgEl(dot, "title", {}, `Turn ${i + 1} · ${show(t[metric])} ${unit} · ${t.verdict || ""} · ${t.note || ""}`);
    });
    svgEl(svg, "text", {x: (left + right) / 2, y: height - 3, "text-anchor": "middle"}, "Turn");
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

    const picker = root.querySelector("#dbsnp-backend");
    data.branches.forEach((br, i) => picker.add(new Option(`${br.branch} (${br.turns} turns)`, i)));
    const withTurns = data.branches.findIndex(br => br.turns > 0 && br.backend !== "snplocs");
    picker.value = String(Math.max(0, withTurns));
    const runtime = root.querySelector("[data-dbsnp-runtime]"), ram = root.querySelector("[data-dbsnp-ram]");
    function paint() {
      const br = data.branches[Number(picker.value)];
      const turns = data.turns.filter(t => t.backend === br.backend && t.kind === "screen");
      chart(runtime, turns, "cold_s", "s", b ? b.cold_s : null);
      chart(ram, turns, "cold_gb", "GB", b ? b.cold_gb : null);
    }
    picker.addEventListener("change", paint);
    new ResizeObserver(() => requestAnimationFrame(paint)).observe(runtime);
    paint();

    const log = root.querySelector("[data-dbsnp-turns] tbody");
    const names = Object.fromEntries(data.branches.map(br => [br.backend, br.branch]));
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
