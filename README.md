# ArtsySelfAssembly

Research code, simulation data, and a static showcase website for a study of
how particles with complex, randomly-heterogeneous patch interactions
self-assemble into crystals, fibers, gels, and other structures.

Paper: *"How Do Particles with Complex Interactions Self-Assemble?"*,
[Phys. Rev. X 14, 041061 (2024)](https://doi.org/10.1103/PhysRevX.14.041061).

## Layout

- **`codes/`** — the original research code: simulation/analysis tools
  (`System.py`, `LatticeTools.py`, `SideClasses.py`) and the lattice geometry
  definitions (`Lattice/`) used by both the original analysis and the site's
  data pipeline below.
- **`data/240531RandomOne/`** — simulation output: precomputed structural
  metrics and category labels for ~9000 simulated assemblies
  (`all_data.pkl`, `labels.pkl`, `labels_unclassified.pkl`). The raw
  per-particle position/orientation data (`ImagesData/`, ~7.5GB) has been
  extracted into the site below and removed to avoid duplicating it.
- **`site/`** — **PatchWork**, a static, no-backend website for browsing and
  exploring the assemblies: a gallery, a category browser, a phase-diagram
  map, and an interactive per-assembly viewer. See `site/README.md` to run it
  locally or regenerate its data, and `site/ARCHITECTURE.md` for a full
  code walkthrough (useful if you want to change how something is computed
  or rendered later).

## Quick start (the website)

```bash
cd site/public
python3 -m http.server 8000
# open http://localhost:8000/
```
