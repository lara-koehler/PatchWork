"""
Regenerates just the thumbnail images from the already-extracted manifest +
binary shards (site/public/data/). Does not touch the raw 7.5GB simulation
pickles at all -- use this when only the thumbnail rendering changes.

Draws the real patchy-hexagon representation (six viridis-colored wedges per
particle, matching codes/LatticeTools.py's plot_2Dcolored_particle) at 2x
supersampling, downsampled for anti-aliasing.

Run: python3 site/build/regen_thumbnails.py
"""
import json
import os

import numpy as np
from PIL import Image, ImageDraw

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
OUT_DIR = os.path.join(ROOT, "site", "public", "data")
SHARD_DIR = os.path.join(OUT_DIR, "shards")
THUMB_DIR = os.path.join(OUT_DIR, "thumbs")

N_PARTICLES = 500
RECORD_SIZE = N_PARTICLES * 2 * 4 + N_PARTICLES * 1 + N_PARTICLES * 4  # 6500

THUMB_W, THUMB_H = 216, 188
SS = 2  # supersampling factor for anti-aliasing

VERTEX_RATIO = 0.5 / np.cos(np.pi / 6)
PARTICLE_SIZE = 0.95
VERTEX_R = VERTEX_RATIO * PARTICLE_SIZE


def vertex_angle(k):
    return (np.pi / 3) * k + np.pi / 6


VX = [np.cos(vertex_angle(k)) for k in range(6)]
VY = [np.sin(vertex_angle(k)) for k in range(6)]

# Exact port of ParticleRepresentation.myColors1[0] = cm.viridis(linspace(0, 0.9, 6))
PATCH_RGB = [
    (0x44, 0x01, 0x54), (0x42, 0x3d, 0x84), (0x2d, 0x6e, 0x8e),
    (0x1e, 0x99, 0x8a), (0x4d, 0xc2, 0x6b), (0xbd, 0xde, 0x26),
]


def face_color_index(orientation, edge_index):
    # Exact port of plot_2Dcolored_particle's color_face formula.
    return (edge_index + 1 - orientation) % 6


def render_thumbnail(xs, ys, orientations, Lx, Ly, path):
    scale = (THUMB_W * SS) / 62.0
    img = Image.new("RGB", (THUMB_W * SS, THUMB_H * SS), (255, 255, 255))
    draw = ImageDraw.Draw(img)
    ox = (THUMB_W * SS - Lx * scale) / 2
    oy = (THUMB_H * SS - Ly * scale) / 2
    for x, y, o in zip(xs, ys, orientations):
        cx = x * scale + ox
        cy = THUMB_H * SS - (y * scale + oy)
        o = int(o)
        for e in range(6):
            x1 = cx + VX[e] * VERTEX_R * scale
            y1 = cy - VY[e] * VERTEX_R * scale
            e2 = (e + 1) % 6
            x2 = cx + VX[e2] * VERTEX_R * scale
            y2 = cy - VY[e2] * VERTEX_R * scale
            draw.polygon([(cx, cy), (x1, y1), (x2, y2)], fill=PATCH_RGB[face_color_index(o, e)])
    img = img.resize((THUMB_W, THUMB_H), Image.LANCZOS)
    img.save(path, "WEBP", quality=85, method=4)


def main():
    with open(os.path.join(OUT_DIR, "manifest.json")) as f:
        manifest = json.load(f)
    assemblies = manifest["assemblies"]

    by_shard = {}
    for a in assemblies:
        by_shard.setdefault(a["shard"], []).append(a)

    total = 0
    for shard_name, items in sorted(by_shard.items()):
        with open(os.path.join(SHARD_DIR, shard_name), "rb") as f:
            buf = f.read()
        for a in items:
            off = a["block"] * RECORD_SIZE
            xy = np.frombuffer(buf, dtype=np.float32, count=1000, offset=off)
            orient = np.frombuffer(buf, dtype=np.uint8, count=500, offset=off + 4000)
            xs, ys = xy[0::2], xy[1::2]
            out_path = os.path.join(THUMB_DIR, os.path.dirname(a["thumb"]).split("/")[-1], os.path.basename(a["thumb"]))
            render_thumbnail(xs, ys, orient, a["Lx"], a["Ly"], out_path)
            total += 1
        print(f"{shard_name}: {len(items)} thumbnails ({total}/{len(assemblies)} total)")

    print(f"Done: {total} thumbnails regenerated")


if __name__ == "__main__":
    main()
