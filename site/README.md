# Patchy Assembly site

A static, no-backend showcase site for the patchy-particle self-assembly simulations
in `data/240531RandomOne/`.

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

That's it — no npm install, no bundler. The only network dependency is loading
three.js from a CDN (unpkg.com) for the 3D viewer on the assembly detail page.

## Regenerating the data

Only needed if the raw simulation data changes. Requires `numpy` and `Pillow`
(no scipy/matplotlib):

```bash
python3 site/build/extract_data.py
```

This reads ~7.5GB of raw pickles from `data/240531RandomOne/ImagesData/` and
writes ~140MB into `site/public/data/` (manifest + binary position shards +
WebP thumbnails). Takes about 3 minutes. The raw pickles themselves are never
read by the website and should not be committed to git.

## Deploying

`site/public/` is a complete static site — push it as-is to GitHub Pages (or any
static host). `site/build/`, and everything outside `site/`, does not need to be
deployed.
