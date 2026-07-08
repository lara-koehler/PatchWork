import { loadManifest, loadAssemblyRecord, pickRandom, assemblyUrl } from "./data.js";
import { CATEGORY_INFO, METRIC_INFO, stickinessLabel, diversityLabel } from "./constants.js";
import { createViewer } from "./viewer.js";

const params = new URLSearchParams(window.location.search);
const id = params.get("id");
const isSurprise = params.get("surprise") === "1";

const breadcrumbs = document.getElementById("breadcrumbs");
const loadingNote = document.getElementById("loading-note");
const layout = document.getElementById("detail-layout");

function gaugeHtml(key, value, range) {
  const info = METRIC_INFO[key];
  if (!info) return "";
  const [lo, hi] = range;
  const t = Math.max(0, Math.min(1, (value - lo) / Math.max(1e-9, hi - lo)));
  return `<div class="gauge-row">
    <div class="gauge-label"><span>${info.label}</span><span>${info.higher}</span></div>
    <div class="gauge-track"><div class="gauge-fill" style="width:${(t * 100).toFixed(0)}%"></div></div>
  </div>`;
}

async function init() {
  const manifest = await loadManifest();
  const entry = manifest.assemblies.find((a) => a.id === id) || pickRandom(manifest.assemblies);

  const info = CATEGORY_INFO[entry.categoryName] || CATEGORY_INFO.unclassified;
  document.title = `${info.label} — PatchWork`;

  breadcrumbs.innerHTML = `<a href="index.html">Gallery</a> &rsaquo; <a href="categories.html?cat=${entry.categoryName}">${info.label}</a> &rsaquo; this assembly`;

  if (isSurprise) {
    breadcrumbs.insertAdjacentHTML(
      "afterend",
      `<div class="surprise-banner">🎲 Surprise! This is one example drawn at random from the dataset — hit "Surprise me" again for another.</div>`
    );
  }

  document.getElementById("cat-heading").textContent = info.label;
  document.getElementById("cat-heading").style.color = info.color;
  document.getElementById("cat-blurb").textContent = info.blurb + (entry.verified ? " (human-verified example)" : "");
  document.getElementById("link-category").href = `categories.html?cat=${entry.categoryName}`;
  document.getElementById("link-phasemap").href = `phase-map.html?highlight=${encodeURIComponent(entry.id)}`;

  document.getElementById("recipe-mu").textContent = `${stickinessLabel(entry.mu)} (${entry.mu})`;
  document.getElementById("recipe-sigma").textContent = `${diversityLabel(entry.sigma)} (${entry.sigma})`;

  const gaugeKeys = ["sizeMax", "sphericity", "porosity", "savRatio"];
  document.getElementById("gauges").innerHTML = gaugeKeys
    .map((k) => gaugeHtml(k, entry.metrics[k], manifest.metricRanges[k]))
    .join("");

  document.getElementById("btn-another").addEventListener("click", () => {
    const pool = manifest.assemblies.filter((a) => a.categoryName === entry.categoryName && a.id !== entry.id);
    window.location.href = assemblyUrl(pickRandom(pool.length ? pool : manifest.assemblies));
  });

  const record = await loadAssemblyRecord(entry);

  loadingNote.hidden = true;
  layout.hidden = false;

  const viewer = createViewer(document.getElementById("viewer-wrap"));
  viewer.loadAssembly(record, entry);

  const explainer = document.getElementById("view-explainer");
  const explainerText = {
    patch: "Each particle shows six colored patches. Every possible contact between two patches has its own level of stickiness, chosen independently and at random.",
    orientation: "Particles are now colored by their rotation state — six possible orientations, one color each. Particles with the same color are sitting the same way.",
    energy: "Particles are now colored by their total bond energy with their neighbors, compared to the assembly's average — blue particles are unusually settled, red ones unusually strained.",
  };
  const modeButtons = document.querySelectorAll(".viewer-controls button[data-mode]");
  explainer.textContent = explainerText.patch;

  modeButtons.forEach((btn) => {
    btn.addEventListener("click", () => {
      modeButtons.forEach((b) => b.classList.remove("active"));
      btn.classList.add("active");
      viewer.setMode(btn.dataset.mode);
      explainer.textContent = explainerText[btn.dataset.mode];
    });
  });
}

init().catch((err) => {
  console.error(err);
  loadingNote.textContent = "Something went wrong loading this assembly. Try reloading the page.";
  loadingNote.style.color = "#e66767";
});
