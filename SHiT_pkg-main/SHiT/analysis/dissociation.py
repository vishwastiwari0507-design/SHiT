"""
mdanalysis.analysis.dissociation
=================================
Classify H₂O, OH⁻, H₃O⁺, H⁺, Si-OH, Si-H species in XYZ trajectories.

CLI
---
    mdanalysis dissociation -d traj.xyz -i input.txt -p para.txt [-o results.txt]

input.txt::

    Dimension:
    a = 0 38.241
    b = 0 38.241
    c = 0 57.32

    Region:
    A = 0 38.241
    B = 0 38.241
    C = 27 57.32

    PBC:
    T T F

para.txt::

    Bond length:
    H-O = 0.8 1.2
    H-Si = 1.5 1.8
    O-Si = 1.4 1.7

    Bond angle:
    H-O-H = 90 115
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Set, Tuple

import numpy as np

from ..core.trajectory import read_xyz
from ..core.parsers import BoxInputParser, parse_range
from ..core.geometry import dist_pbc, angle_pbc, minimum_image


# ---------------------------------------------------------------------------
# Para file reader
# ---------------------------------------------------------------------------

def _parse_para(path: str) -> Tuple[Dict, Dict]:
    """Return (bond_lengths, bond_angles)."""
    bond_lengths: Dict[str, Tuple[float, float]] = {}
    bond_angles: Dict[str, Tuple[float, float]] = {}
    section = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            low = line.lower()
            if "bond length" in low or "bond_length" in low:
                section = "bl"; continue
            if "bond angle" in low or "bond_angle" in low:
                section = "ba"; continue
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                try:
                    lo, hi = parse_range(v.strip())
                except ValueError:
                    continue
                if section == "bl":
                    bond_lengths[k] = (lo, hi)
                elif section == "ba":
                    bond_angles[k] = (lo, hi)
    return bond_lengths, bond_angles


def _bl(bond_lengths: Dict, a: str, b: str) -> Tuple[float, float]:
    for key in (f"{a}-{b}", f"{b}-{a}"):
        if key in bond_lengths:
            return bond_lengths[key]
    raise KeyError(f"No bond length for {a}-{b}")


# ---------------------------------------------------------------------------
# Classifier
# ---------------------------------------------------------------------------

class SpeciesClassifier:
    def __init__(self, bond_lengths, bond_angles, box_dims, eff_pbc):
        self.bl = bond_lengths
        self.ba = bond_angles
        self.box = np.asarray(box_dims, float)
        self.pbc = np.asarray(eff_pbc, bool)

    def classify(self, elements: List[str],
                 positions: np.ndarray) -> Dict[str, int]:
        """Return species counts for one frame."""
        h_idx = [i for i, e in enumerate(elements) if e == "H"]
        o_idx = [i for i, e in enumerate(elements) if e == "O"]
        si_idx = [i for i, e in enumerate(elements) if e == "Si"]

        bl_HO = self.bl.get("H-O", (0.9, 1.1))
        bl_HSi = self.bl.get("H-Si", (1.4, 1.7))
        bl_OSi = self.bl.get("O-Si", (1.4, 1.6))
        ba_HOH = self.ba.get("H-O-H", (100.0, 110.0))

        # Build O→H bond lists
        o_h_bonds: Dict[int, List[Tuple[int, float]]] = {o: [] for o in o_idx}
        for o in o_idx:
            for h in h_idx:
                d = dist_pbc(positions[o], positions[h], self.box, self.pbc)
                if bl_HO[0] <= d <= bl_HO[1]:
                    o_h_bonds[o].append((h, d))

        si_h_bonds: Dict[int, List[Tuple[int, float]]] = {s: [] for s in si_idx}
        si_o_bonds: Dict[int, List[Tuple[int, float]]] = {s: [] for s in si_idx}
        for s in si_idx:
            for h in h_idx:
                d = dist_pbc(positions[s], positions[h], self.box, self.pbc)
                if bl_HSi[0] <= d <= bl_HSi[1]:
                    si_h_bonds[s].append((h, d))
            for o in o_idx:
                d = dist_pbc(positions[s], positions[o], self.box, self.pbc)
                if bl_OSi[0] <= d <= bl_OSi[1]:
                    si_o_bonds[s].append((o, d))

        assigned_H: Set[int] = set()
        assigned_O: Set[int] = set()
        counts = dict(H2O=0, OH=0, H3O=0, H_free=0, Si_OH=0, Si_H=0)

        # 1. H3O+
        for o, bonded in o_h_bonds.items():
            if o in assigned_O:
                continue
            avail = [(h, d) for h, d in bonded if h not in assigned_H]
            if len(avail) >= 3:
                avail.sort(key=lambda x: x[1])
                for h, _ in avail[:3]:
                    assigned_H.add(h)
                assigned_O.add(o)
                counts["H3O"] += 1

        # 2. H2O
        ang_lo, ang_hi = ba_HOH
        for o, bonded in o_h_bonds.items():
            if o in assigned_O:
                continue
            avail = [(h, d) for h, d in bonded if h not in assigned_H]
            avail.sort(key=lambda x: x[1])
            found = False
            for i in range(len(avail)):
                if found: break
                for j in range(i + 1, len(avail)):
                    h1, h2 = avail[i][0], avail[j][0]
                    ang = angle_pbc(positions[h1], positions[o], positions[h2],
                                    self.box, self.pbc)
                    if ang_lo <= ang <= ang_hi:
                        assigned_O.add(o); assigned_H.add(h1); assigned_H.add(h2)
                        counts["H2O"] += 1; found = True; break

        # 3. Si-OH
        o_bonded_to_si: Set[int] = {o for s in si_idx for o, _ in si_o_bonds[s]}
        for o in o_bonded_to_si:
            if o in assigned_O:
                continue
            avail = [(h, d) for h, d in o_h_bonds.get(o, []) if h not in assigned_H]
            if avail:
                avail.sort(key=lambda x: x[1])
                assigned_O.add(o); assigned_H.add(avail[0][0])
                counts["Si_OH"] += 1

        # 4. OH-
        for o, bonded in o_h_bonds.items():
            if o in assigned_O or o in o_bonded_to_si:
                continue
            avail = [(h, d) for h, d in bonded if h not in assigned_H]
            if len(avail) == 1:
                assigned_O.add(o); assigned_H.add(avail[0][0])
                counts["OH"] += 1

        # 5. Si-H
        for s in si_idx:
            for h, _ in si_h_bonds[s]:
                if h not in assigned_H:
                    assigned_H.add(h); counts["Si_H"] += 1

        # 6. Free H+
        counts["H_free"] = sum(1 for h in h_idx if h not in assigned_H)

        return counts


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_dissociation(
    trajectory: str,
    input_file: str,
    para_file: str,
    output_file: str = "water_dissociation.txt",
    verbose: bool = False,
) -> List[Dict]:
    """Run water dissociation analysis.

    Returns list of per-frame dicts with species counts.
    """
    bp = BoxInputParser(input_file).parse()
    bond_lengths, bond_angles = _parse_para(para_file)

    box_dims = bp.box_dims
    pbc = np.array(bp.pbc, bool)
    region_min = np.array(bp.region_min)
    region_max = np.array(bp.region_max)

    # Effective PBC
    eff_pbc = np.array([
        pbc[i] and abs(bp.region_dims[i] - bp.box_dims[i]) < 1e-6
        for i in range(3)
    ])

    classifier = SpeciesClassifier(bond_lengths, bond_angles, box_dims, eff_pbc)
    results = []

    for frame in read_xyz(trajectory):
        # Wrap & filter to region
        coords = frame.coords.copy()
        for i in range(3):
            if pbc[i] and box_dims[i] > 0:
                coords[:, i] = (
                    bp.box_min[i]
                    + (coords[:, i] - bp.box_min[i]) % box_dims[i]
                )
        mask = np.all(
            (coords >= region_min) & (coords <= region_max), axis=1
        )
        elems = [frame.elements[i] for i in range(len(frame.elements)) if mask[i]]
        pos = coords[mask]

        counts = classifier.classify(elems, pos)
        h_h3o = counts["H_free"] + counts["H3O"]
        total_dissoc = h_h3o + counts["OH"]
        row = dict(
            frame=frame.index,
            H_H3O=h_h3o, OH=counts["OH"], Total=total_dissoc,
            H2O=counts["H2O"], Si_OH=counts["Si_OH"], Si_H=counts["Si_H"],
        )
        results.append(row)

        if verbose or frame.index % 100 == 0:
            print(
                f"Frame {frame.index:>6d}: "
                f"H2O={counts['H2O']:>4d}  OH={counts['OH']:>3d}  "
                f"H+/H3O+={h_h3o:>3d}  Si-OH={counts['Si_OH']:>3d}  "
                f"Si-H={counts['Si_H']:>3d}"
            )

    _write_dissociation_results(output_file, results, bp, eff_pbc)
    print(f"\nProcessed {len(results)} frames → {output_file}")
    return results


def _write_dissociation_results(path, results, bp, eff_pbc):
    with open(path, "w") as fh:
        fh.write("# Water Dissociation Analysis Results\n")
        fh.write(f"# Box: X=[{bp.box_min[0]:.2f},{bp.box_max[0]:.2f}] "
                 f"Y=[{bp.box_min[1]:.2f},{bp.box_max[1]:.2f}] "
                 f"Z=[{bp.box_min[2]:.2f},{bp.box_max[2]:.2f}]\n")
        fh.write(f"# Region: X=[{bp.region_min[0]:.2f},{bp.region_max[0]:.2f}] "
                 f"Y=[{bp.region_min[1]:.2f},{bp.region_max[1]:.2f}] "
                 f"Z=[{bp.region_min[2]:.2f},{bp.region_max[2]:.2f}]\n")
        fh.write(f"# EffPBC: {eff_pbc.tolist()}\n")
        fh.write("# Frame   H+/H3O+   OH-   Total   H2O   Si-OH   Si-H\n")
        fh.write("# " + "-" * 60 + "\n")
        for r in results:
            fh.write(
                f"{r['frame']:>7d}   {r['H_H3O']:>7d}   {r['OH']:>5d}   "
                f"{r['Total']:>5d}   {r['H2O']:>5d}   {r['Si_OH']:>5d}   "
                f"{r['Si_H']:>5d}\n"
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis dissociation",
        description="Water dissociation species analysis from XYZ trajectory.",
    )
    parser.add_argument("-d", "--trajectory", required=True)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--parameters", required=True)
    parser.add_argument("-o", "--output", default="water_dissociation.txt")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    run_dissociation(
        args.trajectory, args.input, args.parameters,
        args.output, args.verbose,
    )


if __name__ == "__main__":
    main()
