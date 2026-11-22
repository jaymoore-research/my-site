"use strict";

const cohort = document.getElementById("cohort");
const meta = document.getElementById("meta");
const detail = document.getElementById("detail");
const baselineNote = document.getElementById("baseline");
const charts = [
  [document.getElementById("runtime"), "seconds", "s"],
  [document.getElementById("ram"), "gib", "GB"],
];
const colours = {tested:"#65746e", kept:"#23744d", reference:"#405e80"};
const statuses = {tested:"Tested", kept:"Promoted", reference:"Baseline measurement"};
const targets = {
  seconds:[[1/2, "2×"], [1/5, "5×"], [1/10, "10×"]],
  gib:[[.75, "25% less"], [.5, "50% less"], [.25, "75% less"]],
};
const packageNames = {mss:"MungeSumstats only", bsg:"BSgenome only", both:"Both packages"};
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
  const changed = packageNames[point.change];
  const title = `${label} · ${changed || "Baseline or no change"} · ${point.github ? "Open commit on GitHub" : "Commit is local"}`;
  const dot = element(svg, point.github ? "a" : "g", {
    ...link, class:"point", "aria-label":title, "data-index":index,
  });
  const paint = {
    class:"point-mark", fill:changed ? colour : "white", stroke:colour,
    opacity:point.status === "tested" ? .7 : 1,
  };
  if (point.change === "bsg") {
    const d = radius * .9;
    element(dot, "rect", {...paint, x:x-d, y:y-d, width:2*d, height:2*d});
  } else if (point.change === "both") {
    const d = radius * 1.2;
    element(dot, "polygon", {...paint, points:`${x},${y-d} ${x+d},${y} ${x},${y+d} ${x-d},${y}`});
  } else {
    element(dot, "circle", {...paint, cx:x, cy:y, r:radius});
  }
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

function paint() {
  if (!active) return;
  charts.forEach(args => draw(...args));
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
    const response = await fetch("progress.json?v=stable", {cache:"no-cache"});
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
