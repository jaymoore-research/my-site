"use strict";

const colours = {tested:"#65746e", kept:"#23744d", reference:"#405e80"};
const statuses = {tested:"Tested", kept:"Promoted", reference:"Baseline measurement"};
const targets = {
  seconds:[[1/2, "2×"], [1/5, "5×"], [1/10, "10×"]],
  gib:[[.75, "25% less"], [.5, "50% less"], [.25, "75% less"]],
};
const packageNames = {mss:"MungeSumstats only", bsg:"BSgenome only", both:"Both packages"};
const changeDescriptions = {
  ca626a046eefe39e8b85552b29212883939fc750:"Reuse RefSNP ID tables",
  c7db7ea7c4fbcc18347f6a67bbac85ac00bd3d9b:"Free unused genome-build tables",
  "94b229292e99d973b83c419b19d20edaf8ae1c84":"Read GRCh37 ID tables as needed",
  d45e4ddd0d7e9d46d794f9c8157651dc054ac147:"Stream RefSNP IDs from disk",
  "5b3b1e0b66798e1b8e8a4eb5f9b5e08257d971cd":"Stream RefSNP IDs and speed up block reads",
  e89ecd2da8dd3b0e52f18a3e4ee39054f8d06078:"Read genome files directly",
  b77919ae8ebd4dca4c079fbca45683b3b34b20ef:"Read genome sequences in chunks",
};
const commitRoot = "https://github.com/neurogenomics/";
const numbers = new Intl.NumberFormat("en-GB", {maximumFractionDigits:2});
const dates = new Intl.DateTimeFormat("en-GB", {dateStyle:"medium", timeStyle:"short", timeZone:"UTC"});
const pointPrompt = "Focus or hover over a mark for measurements. Open it to view the commit on GitHub.";

function element(parent, name, attributes, content) {
  const item = document.createElementNS("http://www.w3.org/2000/svg", name);
  for (const [key, value] of Object.entries(attributes)) item.setAttribute(key, value);
  if (content !== undefined) item.textContent = content;
  parent.append(item);
  return item;
}

function showPoint(section, panel, index) {
  const point = panel.points[index];
  section.querySelectorAll(".point").forEach(dot => {
    const selected = Number(dot.dataset.index) === index;
    dot.classList.toggle("selected", selected);
  });
  const dependency = point.dependency ? " · BSgenome " + point.dependency : "";
  const change = packageNames[point.change]
    ? ` · ${packageNames[point.change]} · ${numbers.format(point.savedSeconds)} s saved vs ${point.pairedBase}`
    : "";
  const link = point.github ? "" : " · MungeSumstats commit is local";
  section.querySelector("[data-detail]").textContent = `${point.commit}${dependency} · ${numbers.format(point.seconds)} s · ${numbers.format(point.gib * 1.073741824)} GB · ${statuses[point.status]}${change}${link} · ${dates.format(new Date(point.date))} UTC`;
}

function drawPoint(svg, x, y, point, index, label, select) {
  const colour = colours[point.status];
  const radius = point.status === "kept" ? 4.3 : 3.4;
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

function drawChanges(svg, notes, changes, left, right, top, bottom, compact) {
  notes.replaceChildren();
  notes.hidden = !compact || changes.length === 0;
  const boxes = [];
  changes.forEach(({x, y, description}, i) => {
    if (compact) {
      element(svg, "text", {x:x+5, y:y-6, class:"change-label"}, i+1);
      const item = document.createElement("li");
      item.textContent = description;
      notes.append(item);
      return;
    }
    const label = element(svg, "text", {class:"change-label"}, description);
    const width = label.getComputedTextLength();
    const preferredX = x + 9 + width <= right ? x + 9 : x - 9 - width;
    const labelX = Math.max(left, Math.min(preferredX, right-width));
    const preferredY = Math.max(top+11, Math.min(y-13, bottom-3));
    const overlapping = boxes.filter(box => labelX+width+8 >= box.left && labelX <= box.right+8);
    const candidates = [preferredY, top+11, bottom-3, ...overlapping.flatMap(box => [box.y-18, box.y+18])];
    const labelY = candidates.sort((a,b) => Math.abs(a-preferredY) - Math.abs(b-preferredY)).find(line =>
      line >= top+11 && line <= bottom-3 && overlapping.every(box => Math.abs(line-box.y) >= 17) &&
      Math.abs(line-y) >= 11
    ) ?? preferredY;
    label.setAttribute("x", labelX);
    label.setAttribute("y", labelY);
    boxes.push({left:labelX, right:labelX+width, y:labelY});
    const endX = Math.max(labelX, Math.min(x, labelX+width));
    element(svg, "line", {x1:x, x2:endX, y1:y, y2:labelY+3, stroke:"#8b9b90", "pointer-events":"none"});
  });
}

function draw(svg, notes, metric, unit, panel, select) {
  const points = panel.points;
  const width = svg.clientWidth;
  const height = svg.clientHeight;
  const left = 53, right = width - 56, top = 28, bottom = height - 34;
  if (right <= left) return;
  svg.replaceChildren();
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
  const baseline = actualValue(panel.baseline, metric);
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
  const changes = [];
  points.forEach((point, i) => {
    if (point.status !== "kept" && !(i === 0 && point.status === "reference")) return;
    const next = actualValue(point, metric);
    if (next >= best) return;
    const drop = 1 - next / best;
    const description = changeDescriptions[point.commitSha];
    if (Number.isFinite(best) && drop >= .20 && description) changes.push({x:x(i), y:y(next), description});
    path += best === Infinity ? `M ${x(i)} ${y(next)}` : ` H ${x(i)} V ${y(next)}`;
    best = next;
  });
  if (path) element(svg, "path", {d:path + ` H ${right}`, fill:"none", stroke:colours.kept, "stroke-width":2});

  points.forEach((point, i) => {
    const amount = actualValue(point, metric);
    const relative = metric === "seconds" ? numbers.format(baseline / amount) + "×"
      : numbers.format(100 * (1 - amount / baseline)) + "% RAM reduction";
    const label = `${point.commit}, ${numbers.format(amount)} ${unit}, ${relative}, ${statuses[point.status]}`;
    drawPoint(svg, x(i), y(amount), point, i, label, select);
  });
  drawChanges(svg, notes, changes, left, right, top, bottom, width < 620);
}

function setup(section, panels, initial) {
  const cohort = section.querySelector("select");
  const meta = section.querySelector("[data-meta]");
  if (!panels.length) {
    meta.textContent = "No benchmark data is available.";
    return;
  }
  const baselineNote = section.querySelector("[data-baseline]");
  const detail = section.querySelector("[data-detail]");
  const charts = [
    [section.querySelector("[data-runtime]"), section.querySelector("[data-runtime-notes]"), "seconds", "s"],
    [section.querySelector("[data-ram]"), section.querySelector("[data-ram-notes]"), "gib", "GB"],
  ];
  let panel, selectedIndex = null, resizeFrame;
  const select = index => {
    selectedIndex = index;
    showPoint(section, panel, index);
  };
  function paint() {
    charts.forEach(([svg, notes, metric, unit]) => draw(svg, notes, metric, unit, panel, select));
    if (selectedIndex !== null) select(selectedIndex);
  }
  function switchCohort() {
    panel = panels[Number(cohort.value)];
    selectedIndex = null;
    const {commit, seconds} = panel.baseline;
    baselineNote.textContent = `Target baseline ${commit} · ${numbers.format(seconds)} s · ${numbers.format(actualValue(panel.baseline, "gib"))} GB`;
    const latest = panel.points.map(point => point.date).sort().at(-1);
    meta.textContent = `${panel.points.length} measured builds · Latest measurement ${dates.format(new Date(latest))} UTC`;
    detail.textContent = pointPrompt;
    paint();
  }
  const counts = new Map();
  panels.forEach((choice, i) => {
    const base = `${choice.kind} · ${choice.studies} files · ${choice.parallel ? "side-by-side" : "serial"}`;
    const n = (counts.get(base) || 0) + 1;
    counts.set(base, n);
    cohort.add(new Option(n === 1 ? base : `${base} (panel ${n})`, i));
  });
  cohort.value = String(Math.max(0, panels.findIndex(initial)));
  section.querySelector("[data-cohort-picker]").hidden = panels.length < 2;
  cohort.disabled = panels.length < 2;
  cohort.addEventListener("change", switchCohort);
  switchCohort();
  new ResizeObserver(() => {
    cancelAnimationFrame(resizeFrame);
    resizeFrame = requestAnimationFrame(paint);
  }).observe(charts[0][0]);
}

async function start() {
  const screening = document.getElementById("screening");
  const holdout = document.getElementById("holdout");
  try {
    const response = await fetch("progress.json?v=reviewed", {cache:"no-cache"});
    if (!response.ok) throw new Error("Could not load benchmark data");
    const {panels} = await response.json();
    setup(screening, panels.filter(p => p.kind === "Screening"), p => p.studies === 3 && p.parallel);
    setup(holdout, panels.filter(p => p.kind === "Holdout"), p => p.kind === "Holdout");
  } catch (error) {
    for (const section of [screening, holdout]) {
      section.querySelector("[data-meta]").textContent = "Benchmark data is temporarily unavailable.";
    }
    console.error(error);
  }
}

start();
