# Architecture

How this site is put together, for future-you. Two halves: a **data pipeline**
(Python, runs on your machine, never shipped) and a **website** (plain HTML/CSS/JS,
this is what's deployed). Read `README.md` first for the quick-start; this file
is the deeper "what does each file do and how do I change X" reference.

## The big picture

```
data/240531RandomOne/*.pkl  --[site/build/*.py]-->  site/public/data/  --[fetch() at runtime]-->  site/public/*.html+js
   (raw sim output,                                  (manifest.json +                              (the actual website,
    not shipped)                                       shards/*.bin +                                 what you deploy)
                                                        thumbs/*.webp)
```

The Python pipeline runs **once** (or whenever the underlying data/rendering
logic changes) and writes lean, static files into `site/public/data/`. The
website only ever reads those already-processed files — it never touches the
raw pickles, never does any computation-heavy work, and needs no server.

---

## Part 1: the data pipeline (`site/build/`)

### Where the raw data used to be

`data/240531RandomOne/all_data.pkl` and `labels.pkl` / `labels_unclassified.pkl`
are still there (small, tens of MB). The 7.5GB `ImagesData/` folder (raw
per-particle positions/orientations) and the unused `all_data_frustration.pkl`
were removed after extraction, since everything they could supply is already
baked into `site/public/data/`. **This means `extract_data.py` (the original,
full pipeline) can no longer be re-run as-is** — it needs `ImagesData/` to get
positions/orientations. If you ever need to change something that requires
re-reading the raw simulation output (e.g. add a wholly new dataset, or change
the metrics pulled from `all_data.pkl`), you'd need that folder back.

Two lighter scripts, `regen_thumbnails.py` and `regen_energy.py`, deliberately
do **not** need `ImagesData/` — they read positions/orientations back out of
the already-extracted `site/public/data/shards/*.bin` files instead. Use these
for anything that only touches rendering or the energy calculation.

### `triangular_lattice.py`

A minimal, dependency-light reimplementation of the parts of
`codes/LatticeTools.py::ReadLatticeParticle` needed here. Deliberately doesn't
import `LatticeTools.py` itself, since that pulls in `scipy`/`matplotlib` at
import time for functionality (LEL design, plotting) this pipeline doesn't
need.

- `TriangularLattice`: loads `codes/Lattice/Triangular.json` (the 6 neighbor
  directions, 60° apart) and `Triangular_structures_1particles.txt` (a lookup
  table `(orientation1, orientation2, edge) -> structure id`). This structure
  id is what indexes into a LEL (Local Energy Landscape) vector to get a bond's
  energy.
- `per_particle_bond_energy(...)`: for one assembly's positions + orientations,
  finds every real lattice-neighbor contact and sums the bond energy each
  particle is party to. **This required a periodic-boundary shear correction**:
  the triangular lattice's rows are staggered by half a lattice spacing, so
  wrapping vertically across the box also has to shift x by half the box width
  for every full box-height crossed. A naive independent-axis (rectangular)
  wrap gets this wrong right at the top/bottom seam — verified empirically
  (produces a spurious bimodal bond-count distribution there) before and after
  the fix. If you ever touch this function, re-run that kind of check before
  trusting the output.

### `extract_data.py` (the original full pipeline — currently not runnable, see above)

For each of the 45 `(mu, sigma)` parameter points: loads the corresponding
`ImagesData/*.pkl` file (all 1000 realizations), and for the ~200 of those that
have entries in `all_data.pkl`, pulls out:
- positions, orientations (from `ImagesData`)
- category label — `labels.pkl` if present there (`verified: true`), else
  `labels_unclassified.pkl` (`verified: false`)
- precomputed metrics (`sizeMax`, `porosity`, `sphericity`, etc.) from
  `all_data.pkl`
- a per-particle bond-energy ("asymmetry") array via `triangular_lattice.py`

...and writes three kinds of output (see "Output format" below). This is the
script to extend if you ever add a wholly new dataset/parameter sweep.

### `regen_thumbnails.py`

Rereads positions/orientations straight out of `site/public/data/shards/*.bin`
(no raw data needed) and redraws every thumbnail. Draws the **real**
patchy-hexagon representation — six triangular wedges per particle, colored
from `cm.viridis(linspace(0, 0.9, 6))`, exactly matching
`ParticleRepresentation.plot_2Dcolored_particle` in `codes/LatticeTools.py` —
at 2x supersampling (drawn oversized, then downsampled with Lanczos
resampling) for anti-aliasing, saved as WebP at quality 85. ~90 seconds for all
9000. Run this if you change thumbnail size, colors, or the particle
representation.

### `regen_energy.py`

Rereads positions/orientations out of the shards, recomputes the per-particle
bond energy via `triangular_lattice.py`, and overwrites just the energy
("asymmetry") bytes in each shard record in place. Needs `all_data.pkl` for the
LEL vectors (still present) but not the raw positions data. Run this if you
change the energy/bond calculation.

### Output format (`site/public/data/`)

- **`manifest.json`**: one JSON object loaded once by every page.
  - `categories`: the 8 category names, index-aligned with the integer labels
    used in the original `labels*.pkl` files (`0` = monomer, ... `7` = liquid).
    This is the *data* mapping — don't confuse it with `CATEGORY_NAMES` in
    `constants.js`, which is a separate, purely cosmetic *display order* (see
    below) that's safe to reorder freely.
  - `muValues`, `sigmaValues`: the 5 and 9 distinct parameter values.
  - `metricRanges`: min/max per metric across all 9000 assemblies, used to
    normalize the gauge bars on the detail page.
  - `assemblies`: array of ~9000 entries, each `{id, shard, block, thumb, mu,
    sigma, idx, category, categoryName, verified, Lx, Ly, metrics: {...}}`.
    `shard`+`block` is how the frontend finds this assembly's binary record
    (see below); `thumb` is a relative path under `data/`.
- **`shards/<paramkey>.bin`**: one binary file per `(mu, sigma)` point (45
  total), each holding 200 fixed-size records back to back (one per
  realization index, `block = idx // 5`). Record layout, 6500 bytes each:
  - bytes `0..4000`: 500× `(x, y)` as `float32`, interleaved
  - bytes `4000..4500`: 500× orientation (`uint8`, values 1-6)
  - bytes `4500..6500`: 500× per-particle energy deviation from the assembly's
    mean ("asymmetry", `float32`)
  - All three sub-ranges start on 4-byte boundaries, so the frontend can view
    them directly as typed arrays without copying.
- **`thumbs/<paramkey>/<idx>.webp`**: one ~216×188 thumbnail per assembly.

---

## Part 2: the website (`site/public/`)

Five static HTML pages, no framework, no bundler, no build step for the HTML/
JS/CSS themselves (only the `data/` folder underneath them is generated).
Pages are linked by plain `<a href>` navigation (not a single-page app) with
query-string state (`?cat=fiber`, `?id=...`, `?highlight=...`) for deep-linking.

### `js/constants.js`

Single source of truth for anything that needs to agree across pages:
- `CATEGORY_NAMES`: **display order only** (monomer, oligomer, fiber, gel,
  sponge, polycrystal, crystal, liquid). Safe to reorder — every assembly
  already carries its resolved `categoryName` string from the manifest, so
  nothing here does an integer→name lookup.
- `CATEGORY_INFO`: label, one-line blurb, and color per category. Colors are
  your own established palette from the original plotting code (darkgreen,
  forestgreen, crimson, orange, gold, dodgerblue, royalblue, darkblue) — kept
  as-is, not re-derived.
  - Note: `unclassified` (assemblies with no confident label at all) isn't one
    of your 8 categories — it gets a neutral gray, not a 9th competing color.
- `PATCH_COLORS`: the 6-step viridis sample, `cm.viridis(linspace(0,0.9,6))`,
  precomputed to hex. This is the *only* palette used for actual patch/face
  identity, in both the Python thumbnail renderer and the JS detail viewer —
  kept identical on purpose.
- `faceColorIndex(orientation, edgeIndex)`: exact port of the `color_face`
  formula in `plot_2Dcolored_particle` — `(edgeIndex + 1 - orientation) mod 6`.
  This is *not* the same as a naive "face = edge + orientation" guess; it's a
  faithful port, verified against the source.
- `vertexAngle(k)`, `VERTEX_RATIO`: the hexagon's vertex geometry (vertices sit
  30° off the neighbor-contact directions), matching
  `ParticleRepresentation.vertices_polar` / `cos_vertices` in `LatticeTools.py`.
- `stickinessLabel(mu)`, `diversityLabel(sigma)`: plain-language bucket labels
  for the two phase-map dials.
- `METRIC_INFO`: label/unit/"higher means" text for each metric shown as a
  gauge on the detail page. (`entropy2`/`entropy7` were removed — their exact
  definition isn't recoverable from the available code, so rather than guess
  at a plain-language description they were dropped rather than mislabeled.)

### `js/data.js`

The only module that talks to `data/`.
- `loadManifest()`: fetches `manifest.json` once per page load, caches the
  promise so multiple callers share one fetch.
- `loadAssemblyRecord(entry)`: fetches just one 6500-byte record out of its
  shard via an HTTP `Range` request (falls back to slicing a full response if
  the host/dev-server ignores `Range`, e.g. `python -m http.server`), returns
  `{xy, orientation, asymmetry}` typed-array views, cached per shard+block.
- `pickRandom(assemblies)`, `assemblyUrl(entry)`: small shared helpers.

### `js/nav.js`

Wires up every `[data-surprise-me]` button (present in the nav on every page)
to lazily load the manifest on first click and jump to a random assembly.
Lazy so pages that don't otherwise need the manifest (the science page) don't
pay for a ~4MB fetch just for a button that's rarely clicked.

### Per-page JS

- **`gallery.js`** (`index.html`): loads the manifest, filters by category/
  size-bucket/verified, sorts (default: shuffle), renders a paginated card grid
  (60 at a time, "Load more" button). All client-side, no server round-trips.
- **`category.js`** (`categories.html`): two modes off one file, chosen by the
  `?cat=` query param. No param → hub view, one section per category (colored
  heading + blurb + preview strip + "View all" link). With `?cat=fiber` → a
  full paginated grid for just that category, same sort semantics as the
  gallery.
- **`phase-map.js`** (`phase-map.html`): a `<canvas>` scatter plot, **not** SVG
  or per-point DOM nodes — with ~9000 points that matters for performance.
  - x-axis = patch diversity (`sigma`, 9 values — the higher-cardinality axis,
    given the wider dimension on purpose). y-axis = stickiness (`mu`, 5
    values).
  - Each `(mu, sigma)` grid cell holds up to 200 realizations, so points are
    jittered within their cell (deterministic PRNG, seeded, so layout is
    stable across reloads) — then **clamped** to stay within the plot's
    padding so jitter can never push a point past the frame or behind the
    axis labels.
  - Hit-testing is a linear scan over currently-visible points on
    `mousemove`/`click` (fast enough at this scale; no spatial index needed).
  - `?highlight=<id>` pins one point with a gold ring and pre-opens its
    tooltip — this is how "See on phase map" from a detail page lands you
    somewhere legible instead of a wall of undifferentiated dots.
  - Tooltip flips from right-of-cursor to left-of-cursor near the right edge
    of the canvas so it's never clipped.
- **`assembly.js`** (`assembly.html`): resolves `?id=`, populates the category
  badge/blurb/recipe/gauges from the manifest, fetches the one binary record it
  needs, hands it to `viewer.js`, and wires the three view-mode buttons.
- **`viewer.js`**: the detail-page particle viewer. **Deliberately 2D, not
  3D** — the simulation itself is flat (z is always 0), so there's nothing
  meaningful to rotate out of plane; this is a pan (drag) + zoom
  (scroll/pinch) canvas, not a WebGL/three.js scene. Three color modes,
  switchable without reloading data:
  - **Patch view** (default): six wedges per hexagon via `faceColorIndex` +
    `PATCH_COLORS` — the real per-patch pattern, not a stand-in.
  - **Orientation view**: whole hexagon flat-filled by
    `PATCH_COLORS[orientation - 1]` — makes same-rotation particles easy to
    spot at a glance.
  - **Energy view**: whole hexagon flat-filled on a blue↔gray↔red diverging
    scale by the particle's bond-energy deviation from the assembly's mean
    (the "asymmetry" field from the shard).
  - Draws a solid black frame at the actual simulation box boundary — not
    decorative, particles can bond across it (periodic boundary condition).

### `css/style.css`

One stylesheet, CSS custom properties at the top (`--bg`, `--accent`, etc.)
for the dark theme. No preprocessor, no CSS-in-JS.

---

## Common changes, where to make them

- **Change a category's color or blurb** → `CATEGORY_INFO` in `constants.js`.
  No data regeneration needed (thumbnails don't bake in category color).
- **Reorder how categories are listed anywhere** → `CATEGORY_NAMES` in
  `constants.js`. Purely cosmetic, safe.
- **Change the patch/particle color palette** → `PATCH_COLORS` in
  `constants.js` (affects the live detail-page viewer immediately) **and**
  `PATCH_RGB` in `site/build/regen_thumbnails.py` (then re-run it — thumbnails
  are baked images, not computed live).
- **Change the energy/bond calculation** → `per_particle_bond_energy` in
  `site/build/triangular_lattice.py`, then run `python3
  site/build/regen_energy.py`.
- **Add/change a gauge on the detail page** → add the metric to
  `METRIC_INFO` in `constants.js` (it must already exist in each assembly's
  `metrics` object in `manifest.json`) and to the `gaugeKeys` array in
  `assembly.js`.
- **Add a wholly new metric not in `all_data.pkl` today** → needs the full
  `extract_data.py` pipeline re-run, which needs `ImagesData/` restored.
- **Change thumbnail size/style** → `site/build/regen_thumbnails.py`, then
  re-run (no raw data needed, ~90s for all 9000).
