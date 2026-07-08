"""
Minimal, dependency-light re-implementation of the parts of
codes/LatticeTools.py::ReadLatticeParticle needed to turn raw particle
positions + orientations into bonds and bond energies.

Deliberately does not import codes/LatticeTools.py: that module pulls in
scipy/matplotlib at import time for functionality (LEL design, plotting)
we don't need here. This file only reconstructs the lattice geometry and
the two-particle structure lookup table from the same source files
(codes/Lattice/Triangular.json and Triangular_structures_1particles.txt).
"""
import json
import os

import numpy as np

LATTICE_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "codes", "Lattice")


class TriangularLattice:
    def __init__(self, directory=LATTICE_DIR):
        with open(os.path.join(directory, "Triangular.json")) as f:
            d = json.load(f)

        self.generators = np.array(d["Generators"])  # (2, 3)
        self.neighbors_lattice = np.array(d["Neighbors"])  # (6, 3)
        self.n_neighbors = d["N_Neighbors"]
        self.n_orientations = d["N_Orientations"]
        # orientations[o-1][edge_index] -> face id shown at that edge for orientation o
        self.orientations = np.array(d["Orientations"])

        # Cartesian offset (unit length) from a particle to its neighbor
        # across edge index e (0-indexed here, edge=e+1 in the 1-indexed
        # convention used by the structure table).
        self.neighbors_cartesian = np.array(
            [self._to_cartesian(n) for n in self.neighbors_lattice]
        )[:, :2]

        self.config_to_structure = self._read_structures(
            os.path.join(directory, "Triangular_structures_1particles.txt")
        )
        # N1: number of "single particle" placeholder structures preceding
        # the real two-particle structures in the full 35-entry table.
        # The LEL vectors stored in all_data.pkl are already cut to
        # LEL_full[N1:], so structure id s maps to LEL_cut[s - N1].
        self.N1 = 7

    def _to_cartesian(self, coord):
        out = np.zeros(3)
        for i in range(len(self.generators)):
            out += coord[i] * self.generators[i]
        return out

    @staticmethod
    def _read_structures(path):
        table = {}
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line == "#END":
                    continue
                o1, o2, edge, s = (int(x) for x in line.split())
                table[(o1, o2, edge)] = s
        return table

    def lel_index(self, s):
        return s - self.N1

    def bond_energy_lookup(self):
        """Build a fast (7,7,7) int array: struct_table[o1,o2,edge] -> lel index (or -1)."""
        table = -np.ones((7, 7, 7), dtype=int)
        for (o1, o2, edge), s in self.config_to_structure.items():
            if o1 > 0 and o2 > 0 and edge > 0:
                table[o1, o2, edge] = self.lel_index(s)
        return table


def per_particle_bond_energy(lattice, xs, ys, orientations, box_lx, box_ly, lel_cut, tol=0.2):
    """
    Returns an (N,) array: sum of bond energies (LEL contributions) over all
    lattice-neighbor contacts found for each particle, using periodic
    wrapping that accounts for the triangular lattice's row shear: each row
    is offset by half a lattice spacing from the one below it (generator2 =
    (0.5, 0.866...)), so wrapping vertically across the box also has to shift
    x by half the box width for every full box-height crossed -- a plain
    independent-axis (rectangular) wrap gets this wrong at the top/bottom
    seam (verified empirically: it produces a spurious bimodal bond-count
    distribution right at that edge; this shear-corrected version doesn't).
    """
    n = len(xs)
    dx = xs[:, None] - xs[None, :]
    dy = ys[:, None] - ys[None, :]
    n_wraps_y = np.round(dy / box_ly)
    dy -= box_ly * n_wraps_y
    dx -= (box_lx / 2) * n_wraps_y
    dx -= box_lx * np.round(dx / box_lx)

    lookup = lattice.bond_energy_lookup()
    energy = np.zeros(n)
    bonds = np.zeros(n, dtype=int)

    o = orientations.astype(int)
    for e in range(lattice.n_neighbors):
        ox, oy = lattice.neighbors_cartesian[e]
        # particle j is neighbor of i across edge e+1 if (pos_j - pos_i) ~ (ox,oy)
        # dx[i,j] = x_i - x_j, so we want -dx ~ ox  =>  dx ~ -ox
        match = (np.abs(-dx - ox) < tol) & (np.abs(-dy - oy) < tol)
        np.fill_diagonal(match, False)
        ii, jj = np.where(match)
        if len(ii) == 0:
            continue
        idx = lookup[o[ii], o[jj], e + 1]
        valid = idx >= 0
        ii, idx = ii[valid], idx[valid]
        np.add.at(energy, ii, lel_cut[idx])
        np.add.at(bonds, ii, 1)

    return energy, bonds
