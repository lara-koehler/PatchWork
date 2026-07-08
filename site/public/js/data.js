import { RECORD_SIZE } from "./constants.js";

let manifestPromise = null;

// Loads data/manifest.json once per page and caches the in-flight/resolved
// promise so multiple callers on the same page share one fetch.
export function loadManifest() {
  if (!manifestPromise) {
    manifestPromise = fetch("data/manifest.json").then((r) => {
      if (!r.ok) throw new Error("failed to load manifest.json");
      return r.json();
    });
  }
  return manifestPromise;
}

const shardCache = new Map();

// Fetches just the one 6500-byte record for `entry` out of its shard file.
// Uses a Range request so a static host only transfers the bytes needed;
// falls back to slicing a full response for hosts/dev-servers that ignore
// Range headers (e.g. `python -m http.server`).
export async function loadAssemblyRecord(entry) {
  const cacheKey = entry.shard + ":" + entry.block;
  if (shardCache.has(cacheKey)) return shardCache.get(cacheKey);

  const offset = entry.block * RECORD_SIZE;
  const url = "data/shards/" + entry.shard;
  const resp = await fetch(url, {
    headers: { Range: `bytes=${offset}-${offset + RECORD_SIZE - 1}` },
  });
  let buf = await resp.arrayBuffer();
  if (resp.status !== 206 || buf.byteLength !== RECORD_SIZE) {
    buf = buf.slice(offset, offset + RECORD_SIZE);
  }

  const record = {
    xy: new Float32Array(buf, 0, 1000),
    orientation: new Uint8Array(buf, 4000, 500),
    asymmetry: new Float32Array(buf, 4500, 500),
  };
  shardCache.set(cacheKey, record);
  return record;
}

export function pickRandom(assemblies) {
  return assemblies[Math.floor(Math.random() * assemblies.length)];
}

export function assemblyUrl(entry) {
  return "assembly.html?id=" + encodeURIComponent(entry.id);
}
