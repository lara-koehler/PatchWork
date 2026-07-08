import { loadManifest, assemblyUrl } from "./data.js";
import { CATEGORY_INFO, CATEGORY_NAMES, stickinessLabel, diversityLabel } from "./constants.js";

const wrap = document.getElementById("phase-map-wrap");
const canvas = document.getElementById("pm-canvas");
const ctx = canvas.getContext("2d");
const tooltip = document.getElementById("pm-tooltip");
const catSelect = document.getElementById("pm-category");
const verifiedCheck = document.getElementById("pm-verified");
const countEl = document.getElementById("pm-count");
const legendEl = document.getElementById("pm-legend");

const PAD = { l: 70, r: 30, t: 20, b: 60 };
const DOT_R = 3.2;
const DOT_R_HOVER = 6;
const DOT_R_PINNED = 7;

// x-axis = patch diversity (sigma, 9 values) -- more values, gets the wider
// dimension. y-axis = average stickiness (mu, 5 values).
let manifest, xValues, yValues, points, visible, cellW, cellH, plotW, plotH, cssW, cssH;
let hovered = null;
let pinnedId = null;
const POINT_MARGIN = 10; // keeps jittered/hovered dots inside the frame

function buildPoints(assemblies) {
  const rand = mulberry32(20240607);
  return assemblies.map((a) => ({
    a,
    xIdx: xValues.indexOf(a.sigma),
    yIdx: yValues.indexOf(a.mu),
    jx: rand() * 2 - 1,
    jy: rand() * 2 - 1,
  }));
}

function mulberry32(seed) {
  let t = seed;
  return function () {
    t |= 0; t = (t + 0x6d2b79f5) | 0;
    let x = Math.imul(t ^ (t >>> 15), 1 | t);
    x = (x + Math.imul(x ^ (x >>> 7), 61 | x)) ^ x;
    return ((x ^ (x >>> 14)) >>> 0) / 4294967296;
  };
}

function layout() {
  cssW = wrap.clientWidth - 40;
  cssH = Math.max(420, Math.min(560, cssW * 0.5));
  const dpr = window.devicePixelRatio || 1;
  canvas.width = cssW * dpr;
  canvas.height = cssH * dpr;
  canvas.style.width = cssW + "px";
  canvas.style.height = cssH + "px";
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  plotW = cssW - PAD.l - PAD.r;
  plotH = cssH - PAD.t - PAD.b;
  // Band layout: each parameter value gets a full cell (plotW / N), with its
  // gridline at the cell's center -- so the frame extends half a cell beyond
  // the outermost values on each axis, leaving room for jitter instead of
  // clamping points onto the outermost gridline itself.
  cellW = plotW / xValues.length;
  cellH = plotH / yValues.length;
}

function gridX(i) { return PAD.l + cellW * (i + 0.5); }
function gridY(i) { return PAD.t + cellH * (yValues.length - 1 - i + 0.5); }

function pointPixel(p) {
  const cx = gridX(p.xIdx);
  const cy = gridY(p.yIdx);
  const jitterR = Math.min(cellW, cellH) * 0.4;
  const x = Math.max(PAD.l + POINT_MARGIN, Math.min(PAD.l + plotW - POINT_MARGIN, cx + p.jx * jitterR));
  const y = Math.max(PAD.t + POINT_MARGIN, Math.min(PAD.t + plotH - POINT_MARGIN, cy + p.jy * jitterR));
  return { x, y };
}

function draw() {
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // frame -- the actual plottable area, half a cell beyond the outermost gridlines
  ctx.strokeStyle = "rgba(255,255,255,0.18)";
  ctx.lineWidth = 1;
  ctx.strokeRect(PAD.l, PAD.t, plotW, plotH);

  // gridlines at each parameter value
  ctx.strokeStyle = "rgba(255,255,255,0.08)";
  xValues.forEach((v, i) => {
    const x = gridX(i);
    ctx.beginPath(); ctx.moveTo(x, PAD.t); ctx.lineTo(x, PAD.t + plotH); ctx.stroke();
  });
  yValues.forEach((v, i) => {
    const y = gridY(i);
    ctx.beginPath(); ctx.moveTo(PAD.l, y); ctx.lineTo(PAD.l + plotW, y); ctx.stroke();
  });

  // axis ticks/labels
  ctx.fillStyle = "#9aa2ad";
  ctx.font = "11px -apple-system, sans-serif";
  ctx.textAlign = "center";
  xValues.forEach((v, i) => {
    ctx.fillText(String(v), gridX(i), PAD.t + plotH + 18);
  });
  ctx.textAlign = "right";
  yValues.forEach((v, i) => {
    ctx.fillText(String(v), PAD.l - 10, gridY(i) + 4);
  });

  // axis titles
  ctx.textAlign = "center";
  ctx.font = "13px Georgia, serif";
  ctx.fillStyle = "#e8eaed";
  ctx.fillText("Patch diversity  (uniform → varied)", PAD.l + plotW / 2, PAD.t + plotH + 42);

  ctx.save();
  ctx.translate(16, PAD.t + plotH / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("Average stickiness  (sticky ←→ repulsive)", 0, 0);
  ctx.restore();

  // points
  visible.forEach((p) => {
    const { x, y } = pointPixel(p);
    const info = CATEGORY_INFO[p.a.categoryName] || CATEGORY_INFO.unclassified;
    const isHover = hovered === p;
    const isPinned = p.a.id === pinnedId;
    ctx.beginPath();
    ctx.arc(x, y, isHover || isPinned ? (isPinned ? DOT_R_PINNED : DOT_R_HOVER) : DOT_R, 0, Math.PI * 2);
    ctx.fillStyle = info.color;
    ctx.globalAlpha = p.a.verified ? 0.95 : 0.55;
    ctx.fill();
    if (isHover || isPinned) {
      ctx.globalAlpha = 1;
      ctx.lineWidth = isPinned ? 2.5 : 1.5;
      ctx.strokeStyle = isPinned ? "#ffd54f" : "#ffffff";
      ctx.stroke();
    }
  });
  ctx.globalAlpha = 1;
}

function applyFilter() {
  const cat = catSelect.value;
  const verifiedOnly = verifiedCheck.checked;
  visible = points.filter((p) => {
    if (cat && p.a.categoryName !== cat) return false;
    if (verifiedOnly && !p.a.verified) return false;
    return true;
  });
  countEl.textContent = `${visible.length.toLocaleString()} recipes shown`;
  draw();
}

function findNearest(mx, my) {
  let best = null, bestD = 14 * 14; // px hit radius
  for (const p of visible) {
    const { x, y } = pointPixel(p);
    const d = (x - mx) ** 2 + (y - my) ** 2;
    if (d < bestD) { bestD = d; best = p; }
  }
  return best;
}

function showTooltipFor(p, mx, my) {
  const info = CATEGORY_INFO[p.a.categoryName] || CATEGORY_INFO.unclassified;
  tooltip.style.display = "block";
  tooltip.innerHTML = `<img src="data/${p.a.thumb}" width="90" style="display:block;border-radius:6px;margin-bottom:6px;background:#f7f6f2;" />
    <strong style="color:${info.color}">${info.label}</strong><br/>
    ${stickinessLabel(p.a.mu)} &middot; ${diversityLabel(p.a.sigma)}<br/>
    Largest structure: ${Math.round(p.a.metrics.sizeMax)} / 500`;

  const tw = tooltip.offsetWidth || 120;
  const gap = 24;
  const left = mx + gap + tw > cssW ? mx - gap - tw : mx + gap;
  tooltip.style.left = Math.max(0, left) + "px";
  tooltip.style.top = my + 14 + "px";
}

function onMove(evt) {
  const rect = canvas.getBoundingClientRect();
  const mx = evt.clientX - rect.left;
  const my = evt.clientY - rect.top;
  const hit = findNearest(mx, my);
  if (hit !== hovered) {
    hovered = hit;
    draw();
  }
  if (hit) {
    showTooltipFor(hit, mx, my);
    canvas.style.cursor = "pointer";
  } else {
    tooltip.style.display = "none";
    canvas.style.cursor = "default";
  }
}

function onClick() {
  if (hovered) window.location.href = assemblyUrl(hovered.a);
}

function buildLegend() {
  legendEl.innerHTML = CATEGORY_NAMES.map((name) => {
    const info = CATEGORY_INFO[name];
    return `<span class="cat-chip" style="min-width:auto;cursor:default;">
      <span><span class="dot" style="background:${info.color}"></span>${info.label}</span>
    </span>`;
  }).join("");
}

async function init() {
  manifest = await loadManifest();
  xValues = manifest.sigmaValues;
  yValues = manifest.muValues;
  points = buildPoints(manifest.assemblies);

  const params = new URLSearchParams(window.location.search);
  pinnedId = params.get("highlight");

  CATEGORY_NAMES.concat(["unclassified"]).forEach((name) => {
    const opt = document.createElement("option");
    opt.value = name;
    opt.textContent = CATEGORY_INFO[name].label;
    catSelect.appendChild(opt);
  });

  buildLegend();
  layout();
  applyFilter();

  if (pinnedId) {
    const pinnedPoint = points.find((p) => p.a.id === pinnedId);
    if (pinnedPoint) {
      const { x, y } = pointPixel(pinnedPoint);
      showTooltipFor(pinnedPoint, x, y);
    }
  }

  catSelect.addEventListener("change", applyFilter);
  verifiedCheck.addEventListener("change", applyFilter);
  canvas.addEventListener("mousemove", onMove);
  canvas.addEventListener("mouseleave", () => {
    hovered = null;
    tooltip.style.display = "none";
    draw();
  });
  canvas.addEventListener("click", onClick);
  window.addEventListener("resize", () => { layout(); draw(); });
}

init().catch((err) => {
  console.error(err);
  countEl.textContent = "Couldn't load the phase map data. Try reloading the page.";
});
