"""
Build-time data pipeline: turns the raw simulation pickles in
data/240531RandomOne/ into the lean, static-hostable dataset consumed by
site/public/.

Never run in the browser / not shipped: only its output (site/public/data/)
is deployed. Requires numpy + Pillow locally (no scipy/matplotlib needed).

Outputs (all under site/public/data/):
  manifest.json          one row per assembly (metrics, category, ids)
  shards/<paramkey>.bin   packed positions+orientations+asymmetry, one file
                          per (mu, sigma) parameter point, one fixed-size
                          record per realization index
  thumbs/<paramkey>/<idx>.webp   small dot-density thumbnail per assembly

Run: python3 site/build/extract_data.py
"""
import json
import os
import pickle
import struct
import sys
import time

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(__file__))
from triangular_lattice import TriangularLattice, per_particle_bond_energy

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
DATA_DIR = os.path.join(ROOT, "data", "240531RandomOne")
IMAGES_DIR = os.path.join(DATA_DIR, "ImagesData")
OUT_DIR = os.path.join(ROOT, "site", "public", "data")
SHARD_DIR = os.path.join(OUT_DIR, "shards")
THUMB_DIR = os.path.join(OUT_DIR, "thumbs")

CATEGORY_NAMES = [
    "monomer", "oligomer", "gel", "polycrystal",
    "fiber", "sponge", "crystal", "liquid",
]

N_PARTICLES = 500
RECORD_SIZE = N_PARTICLES * 2 * 4 + N_PARTICLES * 1 + N_PARTICLES * 4  # 6500 bytes

THUMB_W, THUMB_H = 216, 188  # ~ box aspect ratio (60 x 51.96)
THUMB_SCALE = THUMB_W / 62.0  # small margin around the 60-unit box

DOT_COLORS = {
    1: (229, 57, 53), 2: (251, 140, 0), 3: (253, 216, 53),
    4: (67, 160, 71), 5: (30, 136, 229), 6: (142, 36, 170),
}


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def param_key(mu, sigma):
    return "Af_" + str(int(round(mu * 10))).zfill(3) + "_Ani_" + str(int(round(sigma * 10))).zfill(3)


def image_filename(mu, sigma):
    return os.path.join(IMAGES_DIR, "240531RandomOne_" + param_key(mu, sigma) + ".pkl")


def render_thumbnail(xs, ys, orientations, xlim, ylim, path):
    img = Image.new("RGB", (THUMB_W, THUMB_H), (247, 246, 242))
    draw = ImageDraw.Draw(img)
    ox = -xlim[0] * THUMB_SCALE + (THUMB_W - (xlim[1] - xlim[0]) * THUMB_SCALE) / 2
    oy = -ylim[0] * THUMB_SCALE + (THUMB_H - (ylim[1] - ylim[0]) * THUMB_SCALE) / 2
    r = 2.1
    for x, y, o in zip(xs, ys, orientations):
        px = x * THUMB_SCALE + ox
        py = THUMB_H - (y * THUMB_SCALE + oy)
        color = DOT_COLORS.get(int(o), (120, 120, 120))
        draw.ellipse([px - r, py - r, px + r, py + r], fill=color)
    img.save(path, "WEBP", quality=78, method=4)


def pack_record(xs, ys, orientations, asymmetry):
    buf = bytearray(RECORD_SIZE)
    xy = np.empty(2 * N_PARTICLES, dtype=np.float32)
    xy[0::2] = xs
    xy[1::2] = ys
    buf[0:4000] = xy.tobytes()
    buf[4000:4500] = orientations.astype(np.uint8).tobytes()
    buf[4500:6500] = asymmetry.astype(np.float32).tobytes()
    return bytes(buf)


def main():
    t0 = time.time()
    os.makedirs(SHARD_DIR, exist_ok=True)
    os.makedirs(THUMB_DIR, exist_ok=True)

    lattice = TriangularLattice()

    print("loading all_data.pkl / labels ...")
    all_data = load_pickle(os.path.join(DATA_DIR, "all_data.pkl"))
    labels_verified = load_pickle(os.path.join(DATA_DIR, "labels.pkl"))
    labels_auto = load_pickle(os.path.join(DATA_DIR, "labels_unclassified.pkl"))

    by_param = {}
    for key in all_data:
        mu, sigma, idx = key
        by_param.setdefault((mu, sigma), []).append(idx)

    manifest = []
    n_params = len(by_param)
    for pi, ((mu, sigma), idxs) in enumerate(sorted(by_param.items())):
        pk = param_key(mu, sigma)
        img_path = image_filename(mu, sigma)
        print(f"[{pi+1}/{n_params}] {pk}  ({len(idxs)} realizations)  loading {os.path.basename(img_path)} ...")
        images = load_pickle(img_path)

        thumb_subdir = os.path.join(THUMB_DIR, pk)
        os.makedirs(thumb_subdir, exist_ok=True)

        shard_path = os.path.join(SHARD_DIR, pk + ".bin")
        shard_buf = bytearray(200 * RECORD_SIZE)

        for idx in sorted(idxs):
            key = (mu, sigma, idx)
            d = all_data[key]
            snap = images[idx]

            xs, ys, zs = snap["particles_positions"]
            xs = np.asarray(xs, dtype=np.float64)
            ys = np.asarray(ys, dtype=np.float64)
            orient = np.asarray(snap["particles_orientations"], dtype=np.int64)
            xlim, ylim, zlim = snap["limits"]
            Lx, Ly = xlim[1] - xlim[0], ylim[1] - ylim[0]

            energy, bonds = per_particle_bond_energy(lattice, xs, ys, orient, Lx, Ly, d["LEL"])
            asymmetry = energy - energy.mean()

            verified = key in labels_verified
            category = labels_verified[key] if verified else labels_auto.get(key, -1)

            block = idx // 5
            shard_buf[block * RECORD_SIZE:(block + 1) * RECORD_SIZE] = pack_record(xs, ys, orient, asymmetry)

            thumb_rel = f"thumbs/{pk}/{idx}.webp"
            render_thumbnail(xs, ys, orient, xlim, ylim, os.path.join(THUMB_DIR, pk, f"{idx}.webp"))

            manifest.append({
                "id": f"{pk}_{idx}",
                "shard": pk + ".bin",
                "block": block,
                "thumb": thumb_rel,
                "mu": mu, "sigma": sigma, "idx": idx,
                "category": int(category),
                "categoryName": CATEGORY_NAMES[category] if category >= 0 else "unclassified",
                "verified": bool(verified),
                "Lx": round(Lx, 2), "Ly": round(Ly, 2),
                "metrics": {
                    "sizeMax": round(float(d["sizeMax"]), 3),
                    "fullVolume": round(float(d["fullVolume"]), 3),
                    "totalVolume": round(float(d["totalVolume"]), 3),
                    "porosity": round(float(d["porosity"]), 4),
                    "savRatio": round(float(d["savRatio"]), 4),
                    "sphericity": round(float(d["sphericity"]), 4),
                    "holesPerParticle": round(float(d["holesPerParticle"]), 4),
                    "holeSize": round(float(d["holeSize"]), 4),
                    "entropy2": round(float(d["entropy2"]), 4),
                    "entropy7": round(float(d["entropy7"]), 4),
                },
            })

        with open(shard_path, "wb") as f:
            f.write(shard_buf)

        del images
        elapsed = time.time() - t0
        print(f"    done ({elapsed:.0f}s elapsed)")

    with open(os.path.join(OUT_DIR, "manifest.json"), "w") as f:
        json.dump({"categories": CATEGORY_NAMES, "assemblies": manifest}, f)

    print(f"Wrote {len(manifest)} assemblies in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
