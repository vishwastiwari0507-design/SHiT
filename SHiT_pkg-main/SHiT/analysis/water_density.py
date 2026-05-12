"""
mdanalysis.analysis.water_density
===================================
Z-directional water density profile from XYZ trajectory.

CLI
---
    mdanalysis water-density -d traj.xyz -i input.txt -p para.txt

input.txt::

    Dimension:
    a = 38.2
    b = 38.2
    c = 57

    PBC:
    T T T

    Bin_size = 2

para.txt::

    Bond_length:
    O-H = 0.8-1.2

    Bond_angle:
    H-O-H = 90-115

Output
------
``H2O.txt``             – per-frame oxygen positions of each water molecule  
``density_profile.txt`` – binned density vs z for all frames  
``water_count.txt``     – number of water molecules per frame
"""

from __future__ import annotations

import argparse
import os
from typing import Dict, List, Tuple

import numpy as np

from ..core.trajectory import read_xyz
from ..core.parsers import BoxInputParser, parse_range


# ---------------------------------------------------------------------------
# Para parser
# ---------------------------------------------------------------------------

def _parse_para(path: str) -> Tuple[Tuple, Tuple]:
    """Return (oh_range, hoh_range)."""
    oh = (0.8, 1.2)
    hoh = (90.0, 115.0)
    section = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            low = line.lower()
            if "bond_length" in low or "bond length" in low:
                section = "bl"; continue
            if "bond_angle" in low or "bond angle" in low:
                section = "ba"; continue
            if "=" in line:
                k, v = line.split("=", 1)
                k = k.strip()
                try:
                    lo, hi = parse_range(v.strip())
                except ValueError:
                    continue
                if section == "bl" and "o-h" in k.lower():
                    oh = (lo, hi)
                elif section == "ba" and "h-o-h" in k.lower():
                    hoh = (lo, hi)
    return oh, hoh


# ---------------------------------------------------------------------------
# Water finder
# ---------------------------------------------------------------------------

def find_water_molecules(
    elements: List[str], coords: np.ndarray,
    box: np.ndarray, pbc: np.ndarray,
    oh_range: Tuple, hoh_range: Tuple,
) -> List[Dict]:
    """Return list of water dicts {oxygen, h1, h2, angle, dist1, dist2}."""
    o_idx = [i for i, e in enumerate(elements)
             if e.upper() in ("O", "OW")]
    h_idx = [i for i, e in enumerate(elements)
             if e.upper() in ("H", "HW", "H1", "H2")]

    oh_lo, oh_hi = oh_range
    ang_lo, ang_hi = hoh_range
    used_h: set = set()
    waters = []

    for oi in o_idx:
        nearby = []
        for hi in h_idx:
            if hi in used_h:
                continue
            delta = coords[hi] - coords[oi]
            for k in range(3):
                if pbc[k] and box[k] > 0:
                    delta[k] -= box[k] * np.round(delta[k] / box[k])
            d = float(np.linalg.norm(delta))
            if oh_lo <= d <= oh_hi:
                nearby.append((hi, d))

        if len(nearby) < 2:
            continue

        for i in range(len(nearby)):
            if len(waters) > 0 and nearby[i][0] in used_h:
                continue
            for j in range(i + 1, len(nearby)):
                hi_, hj_ = nearby[i][0], nearby[j][0]
                if hi_ in used_h or hj_ in used_h:
                    continue

                v1 = coords[hi_] - coords[oi]
                v2 = coords[hj_] - coords[oi]
                for k in range(3):
                    if pbc[k] and box[k] > 0:
                        v1[k] -= box[k] * np.round(v1[k] / box[k])
                        v2[k] -= box[k] * np.round(v2[k] / box[k])
                n1, n2 = np.linalg.norm(v1), np.linalg.norm(v2)
                if n1 < 1e-10 or n2 < 1e-10:
                    continue
                cos_a = np.clip(np.dot(v1 / n1, v2 / n2), -1.0, 1.0)
                ang = float(np.degrees(np.arccos(cos_a)))
                if ang_lo <= ang <= ang_hi:
                    used_h.add(hi_); used_h.add(hj_)
                    waters.append(dict(
                        oi=oi, z=float(coords[oi, 2]),
                        angle=ang, d1=nearby[i][1], d2=nearby[j][1],
                    ))
                    break

    return waters


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_water_density(
    trajectory: str,
    input_file: str,
    para_file: str,
    output_dir: str = ".",
) -> Tuple[np.ndarray, np.ndarray]:
    """Run water density analysis.

    Returns (bin_centers, density_matrix) where density_matrix has shape
    (n_frames, n_bins).
    """
    os.makedirs(output_dir, exist_ok=True)

    bp = BoxInputParser(input_file).parse()
    oh_range, hoh_range = _parse_para(para_file)
    _bs = bp.extra("bin_size", 1.0)
    bin_size = float(_bs[1] if isinstance(_bs, tuple) else _bs)

    box = bp.box_dims
    pbc = np.array(bp.pbc, bool)

    z_max = float(bp.box_max[2])
    n_bins = max(1, int(np.ceil(z_max / bin_size)))
    bin_edges = np.linspace(0, z_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_vol = float(box[0]) * float(box[1]) * bin_size
    # Angstroms³ → water density normaliser (~0.0334 molecules/Å³ for bulk water)
    norm = bin_vol * 0.0334 if bin_vol > 0 else 1.0

    h2o_lines = []
    density_rows = []

    for frame in read_xyz(trajectory):
        waters = find_water_molecules(
            frame.elements, frame.coords,
            box, pbc, oh_range, hoh_range,
        )

        # Record oxygen positions
        for w in waters:
            h2o_lines.append(
                f"{frame.index} {frame.coords[w['oi'], 0]:.4f} "
                f"{frame.coords[w['oi'], 1]:.4f} {w['z']:.4f} "
                f"{w['angle']:.2f} {w['d1']:.4f} {w['d2']:.4f}\n"
            )

        # Bin counts
        counts = np.zeros(n_bins)
        for w in waters:
            z = w["z"]
            if pbc[2] and box[2] > 0:
                z = z - box[2] * np.floor(z / box[2])
            bi = min(int(np.floor(z / bin_size)), n_bins - 1)
            counts[bi] += 1
        density_rows.append(counts / norm)

        if frame.index % 50 == 0:
            print(f"  Frame {frame.index}: {len(waters)} H₂O molecules")

    density = np.array(density_rows) if density_rows else np.zeros((0, n_bins))

    # Write H2O.txt
    h2o_path = os.path.join(output_dir, "H2O.txt")
    with open(h2o_path, "w") as fh:
        fh.write("# Frame O_x O_y O_z Angle OH1 OH2\n")
        fh.writelines(h2o_lines)

    # Write density_profile.txt
    dp_path = os.path.join(output_dir, "density_profile.txt")
    with open(dp_path, "w") as fh:
        header = "Bin_center\t" + "\t".join(f"F{i+1}" for i in range(len(density_rows)))
        fh.write(header + "\n")
        for bi, z in enumerate(bin_centers):
            row = f"{z:.3f}" + "".join(f"\t{density[fi, bi]:.6f}"
                                        for fi in range(len(density_rows)))
            fh.write(row + "\n")

    # Write water_count.txt
    wc_path = os.path.join(output_dir, "water_count.txt")
    # Count from H2O.txt lines per frame
    from collections import defaultdict
    frame_counts: Dict[int, int] = defaultdict(int)
    for line in h2o_lines:
        frame_counts[int(line.split()[0])] += 1
    with open(wc_path, "w") as fh:
        fh.write("# Frame  N_water\n")
        for fi in sorted(frame_counts):
            fh.write(f"{fi}  {frame_counts[fi]}\n")

    print(f"\nOutput → {output_dir}/")
    print(f"  H2O.txt, density_profile.txt, water_count.txt")
    return bin_centers, density


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis water-density",
        description="Z-directional water density profile from XYZ trajectory.",
    )
    parser.add_argument("-d", "--trajectory", required=True)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--parameters", required=True)
    parser.add_argument("-o", "--output-dir", default=".")
    args = parser.parse_args(argv)

    run_water_density(args.trajectory, args.input, args.parameters, args.output_dir)


if __name__ == "__main__":
    main()
