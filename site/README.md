# PatchWork

A static, no-backend showcase site for the patchy-particle self-assembly simulations
in `data/240531RandomOne/`. See `../README.md` for the whole-repo overview and
`ARCHITECTURE.md` for a full code walkthrough.

## Layout

- `build/` — the (not-deployed) Python data pipeline. Reads the raw pickles in
  `data/240531RandomOne/` (and the lattice geometry in `codes/Lattice/`) and writes
  the lean, static-hostable dataset into `public/data/`.
- `public/` — the actual website. Everything in here is what you'd deploy (e.g. to
  GitHub Pages) — plain HTML/CSS/JS, no build step required to view it.

## Try it yourself

From `site/public/`, start any static file server and open it in a browser:

```bash
cd site/public
python3 -m http.server 8000
# then open http://localhost:8000/ in your browser
```

That's it — no npm install, no bundler, no external dependencies at all. The
assembly detail viewer is a plain 2D `<canvas>` (pan + zoom) — the simulation
is genuinely flat, so there's nothing 3D to render.

## Regenerating the data

`site/build/extract_data.py` is the full pipeline (raw pickles → manifest +
shards + thumbnails), but it needs `data/240531RandomOne/ImagesData/` (the raw
per-particle positions, ~7.5GB), which was deleted after the initial
extraction since everything it could supply is already baked into
`site/public/data/`. It's only needed again if you want to pull in a new
metric or dataset that isn't already extracted.

Two lighter scripts don't need the raw data at all — they read
positions/orientations back out of the already-extracted
`site/public/data/shards/*.bin`:

```bash
python3 site/build/regen_thumbnails.py   # re-render thumbnails (~90s for all 9000)
python3 site/build/regen_energy.py       # recompute the per-particle bond-energy field
```

Requires `numpy` and `Pillow` (no scipy/matplotlib). See `ARCHITECTURE.md` for
what each script actually does.

## Deploying

`site/public/` is a complete static site — push it as-is to GitHub Pages (or any
static host). `site/build/`, and everything outside `site/`, does not need to be
deployed.
