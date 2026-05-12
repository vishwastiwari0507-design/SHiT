"""
mdanalysis.analysis.rdf
=======================
Radial Distribution Function (RDF) from an XYZ trajectory.

CLI
---
    mdanalysis rdf -d traj.xyz -i input.txt [-o output_dir]

input.txt format::

    Dimension:
    a = 0:38.241
    b = 0:38.241
    c = 0:57.32

    Region:
    A = 0:32.8
    B = 0:32.8
    C = 23:30

    PBC:
    T T F

    RDF_Settings:
    r_max = 10.0
    bin_size = 0.05

    RDF:
    O-O
    O-H
    Si-O

Output
------
``<pair>.txt`` per pair, plus ``all.txt`` combining all pairs.
"""

from __future__ import annotations

import argparse
import os
from itertools import combinations_with_replacement
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..core.trajectory import read_xyz
from ..core.geometry import wrap_coords, dist_matrix_pbc
from ..core.parsers import BoxInputParser

DEFAULT_R_MAX = 10.0
DEFAULT_BIN_SIZE = 0.1


# ---------------------------------------------------------------------------
# Input parsing (extends BoxInputParser)
# ---------------------------------------------------------------------------

def _parse_rdf_input(path: str):
    """Return (box_parser, rdf_pairs, r_max, bin_size)."""
    bp = BoxInputParser(path).parse()

    rdf_pairs: List[Tuple[str, str]] = []
    # extra() can return float, tuple, or None — always coerce safely
    def _scalar(v, default):
        if v is None:
            return default
        if isinstance(v, tuple):
            return v[1]   # "10.0" parsed as (0,10) → take hi
        return float(v)

    r_max = _scalar(bp.extra("r_max"), DEFAULT_R_MAX)
    bin_size = _scalar(bp.extra("bin_size"), DEFAULT_BIN_SIZE)

    section = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            low = line.lower()
            if low.startswith("rdf_setting") or low.startswith("setting"):
                section = "settings"; continue
            if low == "rdf:":
                section = "rdf"; continue
            if section == "settings" and "=" in line:
                k, v = line.split("=", 1)
                k = k.strip().lower()
                if k in ("r_max", "rmax"):
                    r_max = float(v)
                elif k in ("bin_size", "binsize", "dr"):
                    bin_size = float(v)
            elif section == "rdf" and "-" in line:
                parts = line.split("-")
                if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                    rdf_pairs.append((parts[0].strip(), parts[1].strip()))

    return bp, rdf_pairs, r_max, bin_size


# ---------------------------------------------------------------------------
# Core calculator
# ---------------------------------------------------------------------------

class RDFCalculator:
    def __init__(self, box_origin, box_dims, region_origin, region_dims,
                 pbc, r_max=DEFAULT_R_MAX, bin_size=DEFAULT_BIN_SIZE):
        self.box_origin = np.asarray(box_origin, float)
        self.box_dims = np.asarray(box_dims, float)
        self.region_origin = np.asarray(region_origin, float)
        self.region_dims = np.asarray(region_dims, float)
        self.pbc = np.asarray(pbc, bool)
        self.r_max = r_max
        self.bin_size = bin_size

        self.n_bins = int(np.ceil(r_max / bin_size))
        self.actual_r_max = self.n_bins * bin_size
        self.bin_edges = np.linspace(0, self.actual_r_max, self.n_bins + 1)
        self.bin_centers = 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])

        # Effective PBC: False if region does not span the full box
        self.eff_pbc = np.array([
            pbc[i] and abs(region_dims[i] - box_dims[i]) < 1e-8
            for i in range(3)
        ])
        self.region_volume = float(np.prod(self.region_dims))
        self._shell_vols = (4 / 3) * np.pi * (
            self.bin_edges[1:] ** 3 - self.bin_edges[:-1] ** 3
        )

    # ------------------------------------------------------------------
    def _filter_region(self, elements, coords):
        mask = np.ones(len(elements), dtype=bool)
        for i in range(3):
            lo = self.region_origin[i]
            hi = lo + self.region_dims[i]
            mask &= (coords[:, i] >= lo) & (coords[:, i] < hi)
        return [elements[j] for j in range(len(elements)) if mask[j]], coords[mask]

    def _get_element_coords(self, elements, coords, elem):
        idx = [i for i, e in enumerate(elements) if e == elem]
        return coords[idx] if idx else np.empty((0, 3))

    # ------------------------------------------------------------------
    def calculate(self, traj_path: str, pairs: List[Tuple[str, str]],
                  verbose=True) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """Return dict: pair_name -> (r_centers, g_r)."""
        accum = {f"{a}-{b}": (np.zeros(self.n_bins), 0, 0, 0)
                 for a, b in pairs}

        n_frames = 0
        for frame in read_xyz(traj_path):
            coords = wrap_coords(frame.coords, self.box_origin,
                                 self.box_dims, self.pbc)
            elems, coords = self._filter_region(frame.elements, coords)
            if not elems:
                continue
            n_frames += 1

            for e1, e2 in pairs:
                name = f"{e1}-{e2}"
                c1 = self._get_element_coords(elems, coords, e1)
                c2 = self._get_element_coords(elems, coords, e2)
                same = e1 == e2
                if len(c1) == 0 or len(c2) == 0:
                    continue

                D = dist_matrix_pbc(c1, c2, self.box_dims, self.eff_pbc)
                if same:
                    flat = D[np.tril_indices(len(c1), k=-1)]
                else:
                    flat = D.ravel()
                valid = flat[(flat > 0) & (flat < self.actual_r_max)]
                hist, _ = np.histogram(valid, bins=self.bin_edges)

                prev = accum[name]
                accum[name] = (prev[0] + hist,
                               prev[1] + len(c1),
                               prev[2] + len(c2),
                               prev[3] + 1)

            if verbose and n_frames % 200 == 0:
                print(f"  {n_frames} frames processed …")

        results = {}
        for (e1, e2), (name) in zip(pairs, accum):
            name = f"{e1}-{e2}"
            hist, tot_n1, tot_n2, nf = accum[name]
            same = e1 == e2
            if nf == 0 or tot_n1 == 0:
                results[name] = (self.bin_centers.copy(), np.zeros(self.n_bins))
                continue
            avg_n1 = tot_n1 / nf
            avg_n2 = tot_n2 / nf
            rho = (avg_n1 if same else avg_n2) / self.region_volume
            g_r = np.zeros(self.n_bins)
            for i in range(self.n_bins):
                ideal = avg_n1 * rho * self._shell_vols[i] * nf
                if same:
                    ideal /= 2
                if ideal > 0:
                    g_r[i] = hist[i] / ideal
            results[name] = (self.bin_centers.copy(), g_r)

        if verbose:
            print(f"  Total frames: {n_frames}")
        return results


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_rdf(
    trajectory: str,
    input_file: str,
    output_dir: str = ".",
    verbose: bool = True,
) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """Run RDF analysis.

    Returns dict: pair_name -> (r_values, g_r).
    """
    os.makedirs(output_dir, exist_ok=True)

    bp, rdf_pairs, r_max, bin_size = _parse_rdf_input(input_file)

    box_origin = np.array(bp.box_min)
    box_dims = bp.box_dims
    region_origin = np.array(bp.region_min)
    region_dims = bp.region_dims
    pbc = np.array(bp.pbc)

    # Auto-detect elements for all-pairs if none specified
    if not rdf_pairs:
        elements: Set[str] = set()
        for frame in read_xyz(trajectory):
            elements.update(frame.elements)
            if frame.index >= 9:
                break
        rdf_pairs = list(combinations_with_replacement(sorted(elements), 2))

    if verbose:
        print(f"Box:    {box_origin} → {box_origin + box_dims}")
        print(f"Region: {region_origin} → {region_origin + region_dims}")
        print(f"PBC:    {pbc}")
        print(f"r_max={r_max}, bin_size={bin_size}")
        print(f"Pairs:  {[f'{a}-{b}' for a,b in rdf_pairs]}")

    calc = RDFCalculator(
        box_origin, box_dims, region_origin, region_dims,
        pbc, r_max=r_max, bin_size=bin_size,
    )
    results = calc.calculate(trajectory, rdf_pairs, verbose=verbose)

    # --- write individual files ---
    for name, (r, g) in results.items():
        out = os.path.join(output_dir, f"{name}.txt")
        with open(out, "w") as fh:
            fh.write(f"# RDF {name}\n# r_max={r_max}, bin_size={bin_size}\n")
            fh.write("# r(Å)\tg(r)\n")
            for ri, gi in zip(r, g):
                fh.write(f"{ri:.6f}\t{gi:.6f}\n")
        if verbose:
            print(f"  Saved {out}")

    # --- write combined file ---
    if results:
        out = os.path.join(output_dir, "all.txt")
        names = list(results)
        r_ref = results[names[0]][0]
        with open(out, "w") as fh:
            fh.write("# Combined RDF\n")
            fh.write("# r(Å)\t" + "\t".join(f"g(r)_{n}" for n in names) + "\n")
            for i in range(len(r_ref)):
                row = f"{r_ref[i]:.6f}" + "".join(
                    f"\t{results[n][1][i]:.6f}" for n in names
                )
                fh.write(row + "\n")
        if verbose:
            print(f"  Saved {out}")

    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis rdf",
        description="Calculate RDF from an XYZ trajectory.",
    )
    parser.add_argument("-d", "--trajectory", required=True)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output-dir", default=".")
    parser.add_argument("-q", "--quiet", action="store_true")
    args = parser.parse_args(argv)

    run_rdf(args.trajectory, args.input, args.output_dir, verbose=not args.quiet)


if __name__ == "__main__":
    main()
