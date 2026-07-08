import { loadManifest, assemblyUrl } from "./data.js";
import { CATEGORY_NAMES, CATEGORY_INFO } from "./constants.js";

const PAGE_SIZE = 60;
const PREVIEW_SIZE = 8;

const hub = document.getElementById("cat-hub");
const detailWrap = document.getElementById("cat-detail");
const detailGrid = document.getElementById("detail-grid");
const loadMoreBtn = document.getElementById("load-more");
const countEl = document.getElementById("result-count");
const sortSelect = document.getElementById("f-sort");
const titleEl = document.getElementById("cat-title");
const subtitleEl = document.getElementById("cat-subtitle");

function featuredSort(a, b) {
  if (a.verified !== b.verified) return a.verified ? -1 : 1;
  return b.metrics.sizeMax - a.metrics.sizeMax;
}

function cardHtml(a) {
  const info = CATEGORY_INFO[a.categoryName] || CATEGORY_INFO.unclassified;
  const star = a.verified ? '<span class="verified-star" title="Human-verified">&#9733;</span>' : "";
  return `<a class="card" href="${assemblyUrl(a)}">
    <div class="thumb-wrap"><img loading="lazy" src="data/${a.thumb}" alt="${info.label} assembly" width="216" height="188" /></div>
    <div class="meta">
      <span class="cat-badge" style="background:${info.color}22;color:${info.color}">${info.label}</span>
      ${star}
    </div>
  </a>`;
}

function renderHub(all) {
  const byCat = {};
  all.forEach((a) => {
    (byCat[a.categoryName] = byCat[a.categoryName] || []).push(a);
  });

  const html = CATEGORY_NAMES.concat(["unclassified"]).map((name) => {
    const items = (byCat[name] || []).slice().sort(featuredSort);
    if (!items.length) return "";
    const info = CATEGORY_INFO[name];
    const preview = items.slice(0, PREVIEW_SIZE);
    return `<section style="margin-bottom:34px;">
      <div style="display:flex;align-items:baseline;gap:10px;margin-bottom:10px;">
        <h2 style="margin:0;color:${info.color};font-size:1.3rem;">${info.label}</h2>
        <span style="color:var(--text-dim);font-size:0.85rem;">${items.length.toLocaleString()} examples</span>
        <a href="categories.html?cat=${name}" style="margin-left:auto;font-size:0.85rem;color:var(--accent);text-decoration:none;">View all &rarr;</a>
      </div>
      <p style="color:var(--text-dim);margin:0 0 12px;max-width:70ch;">${info.blurb}</p>
      <div class="grid">${preview.map(cardHtml).join("")}</div>
    </section>`;
  }).join("");

  hub.innerHTML = html;
}

let filtered = [];
let shown = 0;

function renderMoreDetail() {
  const next = filtered.slice(shown, shown + PAGE_SIZE);
  detailGrid.insertAdjacentHTML("beforeend", next.map(cardHtml).join(""));
  shown += next.length;
  loadMoreBtn.hidden = shown >= filtered.length;
}

function applyDetailSort() {
  const sort = sortSelect.value;
  switch (sort) {
    case "size-desc":
      filtered.sort((a, b) => b.metrics.sizeMax - a.metrics.sizeMax);
      break;
    case "size-asc":
      filtered.sort((a, b) => a.metrics.sizeMax - b.metrics.sizeMax);
      break;
    case "random":
      for (let i = filtered.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [filtered[i], filtered[j]] = [filtered[j], filtered[i]];
      }
      break;
    default:
      filtered.sort(featuredSort);
  }
  shown = 0;
  detailGrid.innerHTML = "";
  renderMoreDetail();
}

function renderDetail(all, catName) {
  hub.hidden = true;
  detailWrap.hidden = false;
  const info = CATEGORY_INFO[catName] || CATEGORY_INFO.unclassified;
  titleEl.textContent = info.label;
  titleEl.style.color = info.color;
  subtitleEl.textContent = info.blurb;

  filtered = all.filter((a) => a.categoryName === catName);
  countEl.textContent = `${filtered.length.toLocaleString()} assemblies`;
  sortSelect.addEventListener("change", applyDetailSort);
  loadMoreBtn.addEventListener("click", renderMoreDetail);
  applyDetailSort();
}

async function init() {
  const manifest = await loadManifest();
  const all = manifest.assemblies;
  const params = new URLSearchParams(window.location.search);
  const cat = params.get("cat");
  if (cat) {
    renderDetail(all, cat);
  } else {
    renderHub(all);
  }
}

init().catch((err) => {
  console.error(err);
  hub.innerHTML = `<p style="color:#e66767;">Couldn't load the category data. Try reloading the page.</p>`;
});
