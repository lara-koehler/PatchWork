import { loadManifest, assemblyUrl } from "./data.js";
import { CATEGORY_INFO } from "./constants.js";

const PAGE_SIZE = 60;

const grid = document.getElementById("gallery-grid");
const loadMoreBtn = document.getElementById("load-more");
const countEl = document.getElementById("result-count");
const catSelect = document.getElementById("f-category");
const sizeSelect = document.getElementById("f-size");
const sortSelect = document.getElementById("f-sort");
const verifiedCheck = document.getElementById("f-verified");

let all = [];
let filtered = [];
let shown = 0;

function sizeBucket(sizeMax) {
  if (sizeMax < 10) return "tiny";
  if (sizeMax < 100) return "small";
  if (sizeMax < 400) return "medium";
  return "large";
}

function applyFiltersAndSort() {
  const cat = catSelect.value;
  const size = sizeSelect.value;
  const verifiedOnly = verifiedCheck.checked;
  const sort = sortSelect.value;

  filtered = all.filter((a) => {
    if (cat && a.categoryName !== cat) return false;
    if (size && sizeBucket(a.metrics.sizeMax) !== size) return false;
    if (verifiedOnly && !a.verified) return false;
    return true;
  });

  switch (sort) {
    case "size-desc":
      filtered.sort((a, b) => b.metrics.sizeMax - a.metrics.sizeMax);
      break;
    case "size-asc":
      filtered.sort((a, b) => a.metrics.sizeMax - b.metrics.sizeMax);
      break;
    case "round-desc":
      filtered.sort((a, b) => b.metrics.sphericity - a.metrics.sphericity);
      break;
    case "random":
      for (let i = filtered.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [filtered[i], filtered[j]] = [filtered[j], filtered[i]];
      }
      break;
    case "featured":
    default:
      filtered.sort((a, b) => {
        if (a.verified !== b.verified) return a.verified ? -1 : 1;
        return b.metrics.sizeMax - a.metrics.sizeMax;
      });
  }

  shown = 0;
  grid.innerHTML = "";
  countEl.textContent = `${filtered.length.toLocaleString()} assemblies`;
  renderMore();
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

function renderMore() {
  const next = filtered.slice(shown, shown + PAGE_SIZE);
  grid.insertAdjacentHTML("beforeend", next.map(cardHtml).join(""));
  shown += next.length;
  loadMoreBtn.hidden = shown >= filtered.length;
}

async function init() {
  const manifest = await loadManifest();
  all = manifest.assemblies;

  const catCounts = {};
  all.forEach((a) => { catCounts[a.categoryName] = (catCounts[a.categoryName] || 0) + 1; });
  Object.keys(CATEGORY_INFO).forEach((name) => {
    if (!catCounts[name]) return;
    const opt = document.createElement("option");
    opt.value = name;
    opt.textContent = `${CATEGORY_INFO[name].label} (${catCounts[name]})`;
    catSelect.appendChild(opt);
  });

  [catSelect, sizeSelect, sortSelect, verifiedCheck].forEach((el) =>
    el.addEventListener("change", applyFiltersAndSort)
  );
  loadMoreBtn.addEventListener("click", renderMore);

  applyFiltersAndSort();
}

init().catch((err) => {
  console.error(err);
  countEl.textContent = "Couldn't load the gallery data. Try reloading the page.";
});
