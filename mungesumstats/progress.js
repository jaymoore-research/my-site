"use strict";

const cohort = document.getElementById("cohort");
const meta = document.getElementById("meta");
const detail = document.getElementById("detail");
const baselineNote = document.getElementById("baseline");
const milestoneLegend = document.getElementById("milestone-legend");
const charts = [
  [document.getElementById("runtime"), "seconds", "s"],
  [document.getElementById("ram"), "gib", "GiB"],
];
const colours = {tested:"#65746e", kept:"#23744d", reference:"#405e80", noisy:"#8d6727"};
const statuses = {tested:"Tested", kept:"Promoted", reference:"Baseline measurement", noisy:"Noisy"};
const targets = {
  seconds:[[1/2, "2×"], [1/5, "5×"], [1/10, "10×"]],
  gib:[[.75, "25% less"], [.5, "50% less"], [.25, "75% less"]],
};
const markerColour = "#765c96";
const numbers = new Intl.NumberFormat("en-GB", {maximumFractionDigits:2});
const dates = new Intl.DateTimeFormat("en-GB", {dateStyle:"medium", timeStyle:"short", timeZone:"UTC"});
let panels = [], active = null, selectedIndex = null, resizeFrame = null, milestone = -1;

function element(parent, name, attributes, content) {
  const item = document.createElementNS("http://www.w3.org/2000/svg", name);
  for (const [key, value] of Object.entries(attributes)) item.setAttribute(key, value);
  if (content !== undefined) item.textContent = content;
  parent.append(item);
  return item;
}

function select(index) {
  selectedIndex = index;
  const point = active.points[index];
  document.querySelectorAll(".point").forEach(dot => {
    const selected = Number(dot.dataset.index) === index;
    dot.classList.toggle("selected", selected);
    dot.setAttribute("aria-pressed", String(selected));
  });
  const dependency = point.dependency ? " · BSgenome " + point.dependency : "";
  const marker = index === milestone ? " · First optimised BSgenome build" : "";
  detail.textContent = `${point.commit}${dependency} · ${numbers.format(point.seconds)} s · ${numbers.format(point.gib)} GiB · ${statuses[point.status]}${marker} · ${dates.format(new Date(point.date))} UTC`;
}

function draw(svg, metric, unit) {
  const points = active.points;
  const width = svg.clientWidth;
  const height = svg.clientHeight;
  const left = 52, right = width - 70, top = milestone >= 0 ? 22 : 10, bottom = height - 34;
  if (right <= left) return;
  svg.replaceChildren();
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  const max = Math.max(active.baseline[metric], ...points.map(point => point[metric]));
  const step = 10 ** Math.floor(Math.log10(max / 4));
  const interval = [1, 2, 2.5, 5, 10].map(v => v * step).find(v => v >= max / 4);
  const limit = interval * Math.ceil(max / interval);
  const x = i => left + (points.length === 1 ? 0.5 : i / (points.length - 1)) * (right - left);
  const y = value => bottom - value / limit * (bottom - top);

  for (let v = 0; v <= limit + interval / 10; v += interval) {
    element(svg, "line", {x1:left, x2:right, y1:y(v), y2:y(v), stroke:"#e1e8e2"});
    element(svg, "text", {x:left-8, y:y(v)+4, "text-anchor":"end"}, numbers.format(v));
  }
  const tickCount = Math.min(points.length, width < 500 ? 4 : 7);
  for (let t = 0; t < tickCount; t++) {
    const i = Math.round(t * (points.length - 1) / Math.max(tickCount - 1, 1));
    element(svg, "text", {x:x(i), y:bottom+16, "text-anchor":"middle"}, i+1);
  }
  element(svg, "text", {x:(left+right)/2, y:height-3, "text-anchor":"middle"}, "Build order");

  for (const [fraction, label] of targets[metric]) {
    const value = active.baseline[metric] * fraction;
    const guide = element(svg, "line", {
      x1:left, x2:right, y1:y(value), y2:y(value), stroke:"#5d7469", "stroke-dasharray":"4 4",
      class:"target-guide", "data-value":value,
    });
    element(guide, "title", {}, `${label}: ${numbers.format(value)} ${unit}`);
    element(svg, "text", {x:right+6, y:y(value)+4}, label);
  }

  if (milestone >= 0) {
    const markerX = x(milestone);
    element(svg, "line", {
      x1:markerX, x2:markerX, y1:top-4, y2:bottom, stroke:markerColour,
      "stroke-dasharray":"2 3", "stroke-width":1.5, class:"bsgenome-marker",
    });
    element(svg, "text", {x:markerX+6, y:top-8, style:`fill:${markerColour}`}, "BSgenome");
  }

  let path = "", best = Infinity;
  points.forEach((point, i) => {
    if (point.status !== "kept" && !(i === 0 && point.status === "reference")) return;
    if (point[metric] >= best) return;
    const next = point[metric];
    path += best === Infinity ? `M ${x(i)} ${y(next)}` : ` H ${x(i)} V ${y(next)}`;
    best = next;
  });
  if (path) element(svg, "path", {d:path + ` H ${right}`, fill:"none", stroke:colours.kept, "stroke-width":2});

  points.forEach((point, i) => {
    const label = `${point.commit}, ${numbers.format(point[metric])} ${unit}, ${statuses[point.status]}`;
    const dot = element(svg, "circle", {
      cx:x(i), cy:y(point[metric]), r:point.status === "kept" ? 4.3 : 3.4,
      fill:point.status === "noisy" ? "white" : colours[point.status], stroke:colours[point.status],
      opacity:point.status === "tested" ? .66 : 1, tabindex:0, role:"button", class:"point",
      "aria-label":label, "aria-pressed":"false", "data-index":i,
    });
    element(dot, "title", {}, label);
    dot.addEventListener("pointerenter", () => select(i));
    dot.addEventListener("click", () => select(i));
    dot.addEventListener("focus", () => select(i));
    dot.addEventListener("keydown", event => {
      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(i); }
    });
  });
}

function paint() {
  if (!active) return;
  charts.forEach(args => draw(...args));
  if (selectedIndex !== null) select(selectedIndex);
}

function switchCohort() {
  active = panels[Number(cohort.value)];
  selectedIndex = null;
  const points = active.points;
  const stock = active.baseline.stockDependency;
  milestone = stock && points.some(point => point.dependency === stock)
    ? points.findIndex(point => point.measured === "C" && point.dependency && point.dependency !== stock)
    : -1;
  milestoneLegend.hidden = milestone < 0;
  const {commit, seconds, gib} = active.baseline;
  baselineNote.textContent = `Target baseline ${commit} · ${numbers.format(seconds)} s · ${numbers.format(gib)} GiB`;
  const latest = points.map(p => p.date).sort().at(-1);
  meta.textContent = `${points.length} measured builds · ${active.studies} GWAS files · Latest measurement ${dates.format(new Date(latest))} UTC`;
  detail.textContent = "Select a dot for the commit and measurements.";
  paint();
}

async function start() {
  try {
    const response = await fetch("progress.json?v=targets", {cache:"no-cache"});
    if (!response.ok) throw new Error("Could not load benchmark data");
    panels = (await response.json()).panels;
    const counts = new Map();
    panels.forEach((panel, i) => {
      const base = `${panel.kind} · ${panel.studies} files · ${panel.parallel ? "side-by-side" : "serial"}`;
      const n = (counts.get(base) || 0) + 1;
      counts.set(base, n);
      cohort.add(new Option(n === 1 ? base : `${base} (panel ${n})`, i));
    });
    if (!panels.length) throw new Error("No benchmark data");
    cohort.disabled = false;
    switchCohort();
  } catch (error) {
    meta.textContent = "Benchmark data is temporarily unavailable.";
    console.error(error);
  }
}

cohort.addEventListener("change", switchCohort);
new ResizeObserver(() => {
  cancelAnimationFrame(resizeFrame);
  resizeFrame = requestAnimationFrame(paint);
}).observe(charts[0][0]);
start();
