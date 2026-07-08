// Shared constants describing the lattice geometry, category taxonomy and
// color language used across the site. Kept small and dependency-free so
// every page can import it directly.

// Display order only (the integer->name resolution already happened at
// build time in site/build/extract_data.py; this just controls UI ordering).
export const CATEGORY_NAMES = [
  "monomer", "oligomer", "fiber", "gel",
  "sponge", "polycrystal", "crystal", "liquid",
];

// Plain-language framing for each structure category (used in cards / detail page).
// Colors are the author's own established category palette (from the original
// research plotting code) — kept as-is rather than re-derived.
export const CATEGORY_INFO = {
  monomer: {
    label: "Monomers",
    blurb: "Particles that never found a partner — still drifting alone.",
    color: "darkgreen",
  },
  oligomer: {
    label: "Oligomers",
    blurb: "Small clumps of a handful of particles, short of a full structure.",
    color: "forestgreen",
  },
  fiber: {
    label: "Fibers",
    blurb: "Long, thin, rope-like strands built from a repeating link.",
    color: "crimson",
  },
  gel: {
    label: "Gels",
    blurb: "A loose, branching network — connected but far from orderly.",
    color: "orange",
  },
  sponge: {
    label: "Sponges",
    blurb: "A porous, riddled-with-holes mesh of particles.",
    color: "gold",
  },
  polycrystal: {
    label: "Polycrystals",
    blurb: "Several well-ordered patches stitched together with mismatched seams.",
    color: "dodgerblue",
  },
  crystal: {
    label: "Crystals",
    blurb: "A single, cleanly repeating pattern all the way through.",
    color: "royalblue",
  },
  liquid: {
    label: "Liquids",
    blurb: "Fully packed together, but with no long-range order — a dense fluid.",
    color: "darkblue",
  },
  unclassified: {
    label: "Unclassified",
    blurb: "Not yet sorted into one of the named categories.",
    color: "#898781",
  },
};

// The 6 patch colors, indexed by "face id" 0-5 (a fixed patch identity on the
// particle, independent of its current rotation). These are exactly
// cm.viridis(linspace(0, 0.9, 6)) from ParticleRepresentation.myColors1 in
// codes/LatticeTools.py — not a re-invented palette.
export const PATCH_COLORS = ["#440154", "#423d84", "#2d6e8e", "#1e998a", "#4dc26b", "#bdde26"];

// Which face color id shows at rendered edge position e (0-5) for a particle
// with rotation state `orientation` (1-6). Exact port of
// ParticleRepresentation.plot_2Dcolored_particle's
// `color_face = (orientations[1-orientation][i_edge]) % n_faces`, where
// Python's negative-index wraparound on the (1-orientation) row makes this
// equivalent to (edgeIndex + 1 - orientation) mod 6.
export function faceColorIndex(orientation, edgeIndex) {
  return (((edgeIndex + 1 - orientation) % 6) + 6) % 6;
}

// Angle (radians) of hexagon vertex k (0-5), matching
// ParticleRepresentation.vertices_polar = contact_centers_polar + half a step.
export function vertexAngle(k) {
  return (Math.PI / 3) * k + Math.PI / 6;
}

// Vertex-to-center radius for a particle of "size" 1, matching
// ratio_size = 0.5 / cos_vertices[0] used throughout LatticeTools.py.
export const VERTEX_RATIO = 0.5 / Math.cos(Math.PI / 6);

// Binary shard layout produced by site/build/extract_data.py
export const N_PARTICLES = 500;
export const RECORD_SIZE = N_PARTICLES * 2 * 4 + N_PARTICLES * 1 + N_PARTICLES * 4; // 6500

export const METRIC_INFO = {
  sizeMax: { label: "Largest structure", unit: "particles", higher: "Bigger" },
  sphericity: { label: "Compactness", unit: "", higher: "Rounder" },
  porosity: { label: "Porosity", unit: "", higher: "More holes" },
  savRatio: { label: "Surface exposure", unit: "", higher: "More exposed" },
  holesPerParticle: { label: "Hole frequency", unit: "", higher: "More holes" },
};

export function stickinessLabel(mu) {
  if (mu <= -3) return "Very sticky";
  if (mu <= -1) return "Sticky";
  if (mu <= 1) return "Neutral";
  if (mu <= 3) return "Weakly repulsive";
  return "Repulsive";
}

export function diversityLabel(sigma) {
  if (sigma <= 1) return "Nearly uniform patches";
  if (sigma <= 5) return "Some patch variety";
  if (sigma <= 9) return "Mixed patch variety";
  if (sigma <= 13) return "High patch variety";
  return "Extreme patch variety";
}
