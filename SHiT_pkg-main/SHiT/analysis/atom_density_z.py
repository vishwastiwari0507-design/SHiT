"""
mdanalysis.analysis.atom_density_z
=====================================
Atom density / distribution along the Z-axis from an XYZ trajectory.

CLI
---
    mdanalysis atom-density -d traj.xyz -i input.txt -p para.txt

input.txt  – simulation box dimensions, PBC flags, and bin size::

    Dimension:
    a = 0:38.241
    b = 0:38.241
    c = 0:57.32

    PBC:
    T T F

    Bin_size = 0.5

para.txt  – which elements to analyse (one per line)::

    Elements:
    Si
    O
    H

Output
------
``<element>_z_data.txt``            – histogram data (z_center, count)
``<element>_z_distribution.png``    – individual plot per element
``ALL_elements_z_distribution.png`` – combined summary plot
"""

from __future__ import annotations

import argparse
import os
from collections import defaultdict
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from ..core.trajectory import read_xyz
from ..core.parsers import BoxInputParser


# ---------------------------------------------------------------------------
# Para file reader
# ---------------------------------------------------------------------------

def _parse_para(path: str) -> Set[str]:
    """Return set of element symbols to analyse (empty set = all elements)."""
    elements: Set[str] = set()
    in_section = False
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("element"):
                in_section = True
                continue
            if in_section:
                if line.endswith(":"):   # next section header
                    break
                elements.add(line)
    return elements


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def collect_z_data(
    traj_path: str,
    elements_filter: Set[str],
    box_dims: np.ndarray,
    pbc: np.ndarray,
) -> Dict[str, List[float]]:
    """Read XYZ trajectory and collect z-coordinates per element symbol."""
    atom_data: Dict[str, List[float]] = defaultdict(list)
    n_frames = 0

    for frame in read_xyz(traj_path):
        for elem, coord in zip(frame.elements, frame.coords):
            if elements_filter and elem not in elements_filter:
                continue
            z = float(coord[2])
            # Wrap z into the primary box when PBC is on in z
            if pbc[2] and box_dims[2] > 0:
                z = z % box_dims[2]
            atom_data[elem].append(z)

        n_frames += 1
        if n_frames % 50 == 0:
            print(f"  Processing frame {n_frames}...")

    print(f"  Total frames: {n_frames}")
    return dict(atom_data)


def compute_histograms(
    atom_data: Dict[str, List[float]],
    z_min: float,
    z_max: float,
    bin_size: float = 0.5,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, np.ndarray]]:
    """Compute histograms for all elements with a shared z-range and bin size."""
    n_bins = max(1, int(np.ceil((z_max - z_min) / bin_size)))
    bin_edges = np.linspace(z_min, z_max, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    hists: Dict[str, np.ndarray] = {}
    for elem, zlist in atom_data.items():
        h, _ = np.histogram(zlist, bins=bin_edges)
        hists[elem] = h

    return bin_centers, bin_edges, hists


def save_data(
    output_dir: str,
    atom_data: Dict[str, List[float]],
    bin_centers: np.ndarray,
    hists: Dict[str, np.ndarray],
):
    """Write per-element histogram data to text files."""
    os.makedirs(output_dir, exist_ok=True)
    for elem, hist in hists.items():
        path = os.path.join(output_dir, f"{elem}_z_data.txt")
        z_vals = np.array(atom_data[elem])
        with open(path, "w") as fh:
            fh.write(f"# Z-distribution for element: {elem}\n")
            fh.write(f"# Total data points : {len(z_vals)}\n")
            fh.write(f"# Mean={z_vals.mean():.4f}  Std={z_vals.std():.4f}  "
                     f"Min={z_vals.min():.4f}  Max={z_vals.max():.4f}\n")
            fh.write("# z_center\tcount\n")
            for z, h in zip(bin_centers, hist):
                fh.write(f"{z:.6f}\t{h}\n")
        print(f"  Saved {path}")


def plot_distributions(
    output_dir: str,
    atom_data: Dict[str, List[float]],
    bin_centers: np.ndarray,
    hists: Dict[str, np.ndarray],
    show: bool = False,
):
    """Create and save z-distribution plots (requires matplotlib)."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available - skipping plots.")
        return

    os.makedirs(output_dir, exist_ok=True)

    # --- Individual plot per element ---
    for elem in sorted(hists):
        hist = hists[elem]
        z_vals = np.array(atom_data[elem])

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.fill_between(bin_centers, hist, alpha=0.25, color="steelblue")
        ax.plot(bin_centers, hist, color="steelblue", lw=1.8, label=elem)
        ax.scatter(bin_centers, hist, s=18, color="navy", zorder=3, alpha=0.7)
        ax.set_xlabel("Z-coordinate (Angstrom)", fontsize=12)
        ax.set_ylabel("Atom count", fontsize=12)
        ax.set_title(f"Z-distribution of {elem}", fontsize=14)
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=11)
        stats = (f"Mean : {z_vals.mean():.3f} A\n"
                 f"Std  : {z_vals.std():.3f} A\n"
                 f"Min  : {z_vals.min():.3f} A\n"
                 f"Max  : {z_vals.max():.3f} A")
        ax.text(0.98, 0.98, stats, transform=ax.transAxes, fontsize=9,
                va="top", ha="right",
                bbox=dict(boxstyle="round", fc="wheat", alpha=0.75))
        plt.tight_layout()
        out = os.path.join(output_dir, f"{elem}_z_distribution.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        if show:
            plt.show()
        plt.close()
        print(f"  Saved {out}")

    # --- Combined summary plot ---
    if len(hists) > 1:
        colors = plt.cm.tab10(np.linspace(0, 1, len(hists)))
        fig, ax = plt.subplots(figsize=(12, 6))
        for elem, c in zip(sorted(hists), colors):
            ax.plot(bin_centers, hists[elem], color=c, lw=2, label=elem)
        ax.set_xlabel("Z-coordinate (Angstrom)", fontsize=12)
        ax.set_ylabel("Atom count", fontsize=12)
        ax.set_title("Combined Z-distributions", fontsize=14)
        ax.grid(True, alpha=0.3, ls="--")
        ax.legend(fontsize=10)
        plt.tight_layout()
        out = os.path.join(output_dir, "ALL_elements_z_distribution.png")
        plt.savefig(out, dpi=150, bbox_inches="tight")
        if show:
            plt.show()
        plt.close()
        print(f"  Saved {out}")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_atom_density_z(
    trajectory: str,
    input_file: str,
    para_file: str,
    output_dir: str = ".",
    show_plots: bool = False,
) -> Dict[str, np.ndarray]:
    """Run atom density-in-z analysis from an XYZ trajectory.

    Parameters
    ----------
    trajectory  : multi-frame .xyz trajectory file
    input_file  : input.txt with box dimensions, PBC flags, and bin_size
    para_file   : para.txt listing element symbols to analyse
    output_dir  : directory for output files and plots
    show_plots  : display plots interactively

    Returns
    -------
    dict of element_symbol -> histogram counts array
    """
    os.makedirs(output_dir, exist_ok=True)

    # Parse input files
    bp = BoxInputParser(input_file).parse()
    box_dims = bp.box_dims
    pbc = np.array(bp.pbc, bool)
    _bs = bp.extra("bin_size", 0.5)
    bin_size = float(_bs[1] if isinstance(_bs, tuple) else _bs)

    elements_filter = _parse_para(para_file)

    z_min = float(bp.box_min[2])
    z_max = float(bp.box_max[2])

    print(f"Reading {trajectory}...")
    print(f"  Box z-range : [{z_min:.3f}, {z_max:.3f}] Angstrom")
    print(f"  Bin size    : {bin_size} Angstrom")
    print(f"  PBC         : {pbc.tolist()}")
    print(f"  Elements    : {sorted(elements_filter) if elements_filter else 'all'}")

    atom_data = collect_z_data(trajectory, elements_filter, box_dims, pbc)

    if not atom_data:
        print("No matching atoms found. Check element names in para.txt.")
        return {}

    for elem in sorted(atom_data):
        print(f"  {elem}: {len(atom_data[elem])} data points")

    bin_centers, bin_edges, hists = compute_histograms(
        atom_data, z_min, z_max, bin_size
    )
    save_data(output_dir, atom_data, bin_centers, hists)
    plot_distributions(output_dir, atom_data, bin_centers, hists, show=show_plots)

    print(f"\nDone -> {output_dir}/")
    return hists


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis atom-density",
        description="Atom density along Z from an XYZ trajectory.",
    )
    parser.add_argument("-d", "--trajectory", required=True,
                        help="Multi-frame .xyz trajectory file")
    parser.add_argument("-i", "--input", required=True,
                        help="input.txt with box dimensions, PBC, and bin_size")
    parser.add_argument("-p", "--parameters", required=True,
                        help="para.txt listing element symbols to analyse")
    parser.add_argument("-o", "--output-dir", default="z_distribution_plots",
                        help="Output directory (default: z_distribution_plots)")
    parser.add_argument("--show", action="store_true",
                        help="Display plots interactively")
    args = parser.parse_args(argv)

    run_atom_density_z(
        args.trajectory,
        args.input,
        args.parameters,
        output_dir=args.output_dir,
        show_plots=args.show,
    )


if __name__ == "__main__":
    main()
