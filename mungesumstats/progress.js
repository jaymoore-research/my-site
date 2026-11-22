"use strict";

const cohort = document.getElementById("cohort");
const meta = document.getElementById("meta");
const detail = document.getElementById("detail");
const packageFilter = document.getElementById("package");
const savingsMeta = document.getElementById("savings-meta");
const savingsSvg = document.getElementById("savings");
const baselineNote = document.getElementById("baseline");
const charts = [
  [document.getElementById("runtime"), "seconds", "s"],
  [document.getElementById("ram"), "gib", "GB"],
];
const colours = {tested:"#65746e", kept:"#23744d", reference:"#405e80", noisy:"#8d6727"};
const statuses = {tested:"Tested", kept:"Promoted", reference:"Baseline measurement", noisy:"Noisy"};
const targets = {
  seconds:[[1/2, "2×"], [1/5, "5×"], [1/10, "10×"]],
  gib:[[.75, "25% less"], [.5, "50% less"], [.25, "75% less"]],
};
const packageNames = {mss:"MungeSumstats only", bsg:"BSgenome only", both:"Both packages"};
const packageColours = {mss:"#36658e", bsg:"#765c96", both:"#23744d"};
const commitRoot = "https://github.com/neurogenomics/";
const numbers = new Intl.NumberFormat("en-GB", {maximumFractionDigits:2});
const dates = new Intl.DateTimeFormat("en-GB", {dateStyle:"medium", timeStyle:"short", timeZone:"UTC"});
let panels = [], active = null, selectedIndex = null, resizeFrame = null;

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

  });
  const dependency = point.dependency ? " · BSgenome " + point.dependency : "";
  const change = packageNames[point.change]
    ? ` · ${packageNames[point.change]} · ${numbers.format(point.savedSeconds)} s saved vs ${point.pairedBase}`
    : "";
  const link = point.github ? "" : " · MungeSumstats commit is local";
  detail.textContent = `${point.commit}${dependency} · ${numbers.format(point.seconds)} s · ${numbers.format(point.gib * 1.073741824)} GB · ${statuses[point.status]}${change}${link} · ${dates.format(new Date(point.date))} UTC`;
}

function drawPoint(svg, x, y, point, index, label, colour, radius) {
  const link = point.github
    ? {href:`${commitRoot}mungesumstats-autoresearch/commit/${point.commitSha}`, target:"_blank", rel:"noopener noreferrer"}
    : {role:"button", tabindex:0};
  const title = label + (point.github ? " · Open commit on GitHub" : " · Commit is local");
  const dot = element(svg, point.github ? "a" : "g", {
    ...link, class:"point", "aria-label":title, "data-index":index,
  });
  element(dot, "circle", {
    cx:x, cy:y, r:radius, fill:point.status === "noisy" ? "white" : colour, stroke:colour,
    opacity:point.status === "tested" ? .7 : 1,
  });
  element(dot, "title", {}, title);
  dot.addEventListener("pointerenter", () => select(index));
  dot.addEventListener("focus", () => select(index));
  if (!point.github) {
    dot.addEventListener("click", () => select(index));
    dot.addEventListener("keydown", event => {
      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(index); }
    });
  }
}

function actualValue(point, metric) {
  return point[metric] * (metric === "gib" ? 1.073741824 : 1);
}

function draw(svg, metric, unit) {
  const points = active.points;
  const width = svg.clientWidth;
  const height = svg.clientHeight;
  const left = 53, right = width - 56, top = 28, bottom = height - 34;
  if (right <= left) return;
  svg.replaceChildren();
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  const baseline = actualValue(active.baseline, metric);
  const max = Math.max(baseline, ...points.map(point => actualValue(point, metric)));
  const step = 10 ** Math.floor(Math.log10(max / 4));
  const interval = [1, 2, 2.5, 5, 10].map(v => v * step).find(v => v >= max / 4);
  const limit = interval * Math.ceil(max / interval);
  const x = i => left + (points.length === 1 ? 0.5 : i / (points.length - 1)) * (right - left);
  const y = value => bottom - value / limit * (bottom - top);

  element(svg, "text", {x:1, y:12}, metric === "seconds" ? "Speedup" : "RAM saved");
  element(svg, "text", {x:right+7, y:12}, unit);
  for (let v = 0; v <= limit + interval / 10; v += interval) {
    element(svg, "line", {x1:left, x2:right, y1:y(v), y2:y(v), stroke:"#e1e8e2"});
    element(svg, "text", {x:right+7, y:y(v)+4}, numbers.format(v));
  }
  const tickCount = Math.min(points.length, width < 500 ? 4 : 7);
  for (let t = 0; t < tickCount; t++) {
    const i = Math.round(t * (points.length - 1) / Math.max(tickCount - 1, 1));
    element(svg, "text", {x:x(i), y:bottom+16, "text-anchor":"middle"}, i+1);
  }
  element(svg, "text", {x:(left+right)/2, y:height-3, "text-anchor":"middle"}, "Build order");

  for (const [fraction, label] of targets[metric]) {
    const value = baseline * fraction;
    const guide = element(svg, "line", {
      x1:left, x2:right, y1:y(value), y2:y(value), stroke:"#5d7469", "stroke-dasharray":"4 4",
      class:"target-guide", "data-value":value,
    });
    element(guide, "title", {}, `${label}: ${numbers.format(value)} ${unit}`);
    element(svg, "text", {x:left-6, y:y(value)+4, "text-anchor":"end"}, label.replace(" less", ""));
  }

  let path = "", best = Infinity;
  points.forEach((point, i) => {
    if (point.status !== "kept" && !(i === 0 && point.status === "reference")) return;
    const next = actualValue(point, metric);
    if (next >= best) return;
    path += best === Infinity ? `M ${x(i)} ${y(next)}` : ` H ${x(i)} V ${y(next)}`;
    best = next;
  });
  if (path) element(svg, "path", {d:path + ` H ${right}`, fill:"none", stroke:colours.kept, "stroke-width":2});

  points.forEach((point, i) => {
    const amount = actualValue(point, metric);
    const relative = metric === "seconds" ? numbers.format(baseline / amount) + "×"
      : numbers.format(100 * (1 - amount / baseline)) + "% RAM reduction";
    const label = `${point.commit}, ${numbers.format(amount)} ${unit}, ${relative}, ${statuses[point.status]}`;
    drawPoint(svg, x(i), y(amount), point, i, label, colours[point.status], point.status === "kept" ? 4.3 : 3.4);
  });
}

function drawSavings() {
  const points = active.points;
  const included = points.map((point, index) => [point, index]).filter(([point]) =>
    packageNames[point.change] && (packageFilter.value === "all" || point.change === packageFilter.value));
  const width = savingsSvg.clientWidth, height = savingsSvg.clientHeight;
  const left = 52, right = width - 56, top = 15, bottom = height - 34;
  savingsSvg.replaceChildren();
  savingsSvg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  const counts = Object.keys(packageNames).map(key => [
    key, included.filter(([point]) => point.change === key).length,
  ]).filter(([, count]) => count);
  savingsMeta.textContent = counts.map(([key,count]) => `${packageNames[key]}: ${count}`).join(" · ") || "No measured comparisons in this group.";
  if (!included.length || right <= left) return;

  const values = included.map(([point]) => point.savedSeconds);
  const low = Math.min(0, ...values), high = Math.max(0, ...values);
  const span = high - low || 1;
  const step = 10 ** Math.floor(Math.log10(span / 4));
  const interval = [1, 2, 2.5, 5, 10].map(v => v * step).find(v => v >= span / 4);
  const min = interval * Math.floor(low / interval);
  const max = interval * (high ? Math.ceil(high / interval) : 1);
  const x = i => left + i / Math.max(points.length - 1, 1) * (right - left);
  const y = value => bottom - (value - min) / (max - min) * (bottom - top);
  for (let v = min; v <= max + interval/10; v += interval) {
    element(savingsSvg, "line", {x1:left, x2:right, y1:y(v), y2:y(v), stroke:v === 0 ? "#708476" : "#e1e8e2"});
    element(savingsSvg, "text", {x:right+7, y:y(v)+4}, numbers.format(v));
  }
  element(savingsSvg, "text", {x:(left+right)/2, y:height-3, "text-anchor":"middle"}, "Build order");
  included.forEach(([point, index]) => {
    const label = `${point.commit}, ${packageNames[point.change]}, ${numbers.format(point.savedSeconds)} seconds saved`;
    drawPoint(savingsSvg, x(index), y(point.savedSeconds), point, index, label, packageColours[point.change], 4.2);
  });
}

function paint() {
  if (!active) return;
  charts.forEach(args => draw(...args));
  drawSavings();
  if (selectedIndex !== null) select(selectedIndex);
}

function switchCohort() {
  active = panels[Number(cohort.value)];
  selectedIndex = null;
  const points = active.points;
  const {commit, seconds, gib} = active.baseline;
  baselineNote.textContent = `Target baseline ${commit} · ${numbers.format(seconds)} s · ${numbers.format(gib * 1.073741824)} GB`;
  const latest = points.map(p => p.date).sort().at(-1);
  meta.textContent = `${points.length} measured builds · ${active.studies} GWAS files · Latest measurement ${dates.format(new Date(latest))} UTC`;
  detail.textContent = "Focus or hover over a dot for measurements. Open it to view the commit on GitHub.";
  paint();
}

async function start() {
  try {
    const response = await fetch("progress.json?v=commits", {cache:"no-cache"});
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
packageFilter.addEventListener("change", paint);
new ResizeObserver(() => {
  cancelAnimationFrame(resizeFrame);
  resizeFrame = requestAnimationFrame(paint);
}).observe(charts[0][0]);
start();
