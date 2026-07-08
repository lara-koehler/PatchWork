import { loadManifest } from "./data.js";
import { pickRandom, assemblyUrl } from "./data.js";

// Wires up every element with [data-surprise-me] to jump to a random
// assembly. Manifest is only fetched lazily, on first click, so pages that
// don't otherwise need it (science page, detail page) don't pay for it
// on load.
function initSurpriseMe() {
  const buttons = document.querySelectorAll("[data-surprise-me]");
  buttons.forEach((btn) => {
    btn.addEventListener("click", async () => {
      btn.disabled = true;
      btn.textContent = "Shuffling…";
      try {
        const manifest = await loadManifest();
        const entry = pickRandom(manifest.assemblies);
        window.location.href = assemblyUrl(entry) + "&surprise=1";
      } catch (e) {
        btn.disabled = false;
        btn.textContent = "Surprise me";
        console.error(e);
      }
    });
  });
}

document.addEventListener("DOMContentLoaded", initSurpriseMe);
