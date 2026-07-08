"""
Recomputes the per-particle bond-energy ("asymmetry") field in-place inside
the already-extracted binary shards, using the corrected shear-aware
periodic-boundary handling in triangular_lattice.py. Reads positions and
orientations back out of the shards themselves (not the raw 7.5GB
ImagesData pickles, which aren't needed here) and only needs all_data.pkl
for the LEL vectors.

Run: python3 site/build/regen_energy.py
"""
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from triangular_lattice import TriangularLattice, per_particle_bond_energy

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DATA_DIR = os.path.join(ROOT, "data", "240531RandomOne")
OUT_DIR = os.path.join(ROOT, "site", "public", "data")
SHARD_DIR = os.path.join(OUT_DIR, "shards")

N_PARTICLES = 500
RECORD_SIZE = N_PARTICLES * 2 * 4 + N_PARTICLES * 1 + N_PARTICLES * 4  # 6500


def main():
    lattice = TriangularLattice()

    with open(os.path.join(DATA_DIR, "all_data.pkl"), "rb") as f:
        all_data = pickle.load(f)
    with open(os.path.join(OUT_DIR, "manifest.json")) as f:
        manifest = json.load(f)
    assemblies = manifest["assemblies"]

    by_shard = {}
    for a in assemblies:
        by_shard.setdefault(a["shard"], []).append(a)

    total = 0
    for shard_name, items in sorted(by_shard.items()):
        path = os.path.join(SHARD_DIR, shard_name)
        with open(path, "rb") as f:
            buf = bytearray(f.read())

        for a in items:
            key = (a["mu"], a["sigma"], a["idx"])
            lel = all_data[key]["LEL"]

            off = a["block"] * RECORD_SIZE
            xy = np.frombuffer(buf, dtype=np.float32, count=1000, offset=off).copy()
            orient = np.frombuffer(buf, dtype=np.uint8, count=500, offset=off + 4000).copy()
            xs, ys = xy[0::2], xy[1::2]

            energy, _ = per_particle_bond_energy(lattice, xs, ys, orient, a["Lx"], a["Ly"], lel)
            asymmetry = (energy - energy.mean()).astype(np.float32)
            buf[off + 4500:off + 6500] = asymmetry.tobytes()
            total += 1

        with open(path, "wb") as f:
            f.write(buf)
        print(f"{shard_name}: {len(items)} assemblies re-computed ({total}/{len(assemblies)} total)")

    print(f"Done: {total} assemblies' energy data corrected")


if __name__ == "__main__":
    main()
