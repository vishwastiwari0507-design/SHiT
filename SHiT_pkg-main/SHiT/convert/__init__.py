"""
mdanalysis.convert
==================
Trajectory format conversion and frame extraction utilities.

CLI commands
------------
    mdanalysis lammps-to-xyz  input.lammpstrj output.xyz
    mdanalysis extract-frames input.xyz output.xyz [--frames 0-5,10:50:5,100:500:50]
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, Set


# ---------------------------------------------------------------------------
# LAMMPS → XYZ
# ---------------------------------------------------------------------------

def lammps_to_xyz(trj_path: str, xyz_path: str) -> int:
    """Convert a LAMMPS dump to multi-frame XYZ.

    Requires the dump to have columns: id element x y z
    Returns the number of frames written.
    """
    n_frames = 0
    with open(trj_path) as fin, open(xyz_path, "w") as fout:
        while True:
            line = fin.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = fin.readline().strip()
            fin.readline()                      # ITEM: NUMBER OF ATOMS
            natoms = int(fin.readline().strip())
            fin.readline()                      # ITEM: BOX BOUNDS
            for _ in range(3):
                fin.readline()

            header = fin.readline().strip()
            if not header.startswith("ITEM: ATOMS"):
                raise ValueError("Unexpected LAMMPS format (no ATOMS header)")

            cols = header.split()[2:]
            try:
                ci = {c: i for i, c in enumerate(cols)}
                id_i = ci["id"]; el_i = ci["element"]
                x_i = ci["x"]; y_i = ci["y"]; z_i = ci["z"]
            except KeyError as exc:
                raise ValueError(f"LAMMPS dump missing column: {exc}")

            atoms = []
            for _ in range(natoms):
                parts = fin.readline().split()
                atoms.append((
                    int(parts[id_i]),
                    parts[el_i],
                    float(parts[x_i]),
                    float(parts[y_i]),
                    float(parts[z_i]),
                ))
            atoms.sort(key=lambda a: a[0])

            fout.write(f"{natoms}\nTimestep {timestep}\n")
            for _, el, x, y, z in atoms:
                fout.write(f"{el} {x:.6f} {y:.6f} {z:.6f}\n")
            n_frames += 1

    print(f"Converted {n_frames} frames → {xyz_path}")
    return n_frames


# ---------------------------------------------------------------------------
# Frame extraction
# ---------------------------------------------------------------------------

def parse_frame_spec(spec: str) -> Set[int]:
    """Parse a comma-separated frame specification into a set of indices.

    Supported formats (can be mixed with commas):
    - ``5``         single frame
    - ``0-5``       inclusive range step 1
    - ``10:50:5``   start:stop:step (like Python range, but inclusive stop)
    - ``100:500:50``

    Examples::

        "0-5,10:50:5,100:500:50,1000:5000:500"
    """
    frames: Set[int] = set()
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            tokens = part.split(":")
            start, stop, step = int(tokens[0]), int(tokens[1]), int(tokens[2]) if len(tokens) > 2 else 1
            frames.update(range(start, stop + 1, step))
        elif "-" in part and not part.startswith("-"):
            a, b = part.split("-", 1)
            frames.update(range(int(a), int(b) + 1))
        else:
            frames.add(int(part))
    return frames


def default_frame_list() -> Set[int]:
    """Return the default frame set used in the original scripts."""
    frames: Set[int] = set()
    frames.update(range(0, 51, 5))
    frames.update(range(100, 501, 50))
    frames.update(range(1000, 5001, 500))
    frames.update(range(5000, 25001, 2500))
    return frames


def extract_frames(
    input_path: str,
    output_path: str,
    frames_to_keep: Optional[Set[int]] = None,
    verbose: bool = False,
) -> int:
    """Extract specific frames from an XYZ trajectory.

    Parameters
    ----------
    input_path    : source XYZ file
    output_path   : destination XYZ file
    frames_to_keep: set of 0-based frame indices; None → use default set
    verbose       : print per-frame progress

    Returns number of frames written.
    """
    if frames_to_keep is None:
        frames_to_keep = default_frame_list()

    max_needed = max(frames_to_keep) if frames_to_keep else 0
    extracted = 0
    current = 0

    with open(input_path) as fin, open(output_path, "w") as fout:
        while True:
            line = fin.readline()
            if not line:
                break
            try:
                n_atoms = int(line.strip())
            except ValueError:
                continue

            comment = fin.readline()
            atom_lines = [fin.readline() for _ in range(n_atoms)]

            if current in frames_to_keep:
                fout.write(f"{n_atoms}\n")
                comment_text = comment.strip()
                fout.write(
                    (comment_text + f" | Original frame: {current}" if comment_text
                     else f"Frame {current}") + "\n"
                )
                fout.writelines(atom_lines)
                extracted += 1
                if verbose:
                    print(f"  Extracted frame {current} ({extracted})")

            current += 1

            if current > max_needed and extracted == len(frames_to_keep):
                print("  All requested frames extracted.")
                break

    print(f"\nExtracted {extracted} frames from {current} total → {output_path}")
    return extracted


def count_xyz_frames_quick(path: str) -> int:
    """Count frames without loading coordinates."""
    from ..core.trajectory import count_xyz_frames
    return count_xyz_frames(path)


# ---------------------------------------------------------------------------
# CLI helpers (called by top-level dispatcher)
# ---------------------------------------------------------------------------

def _lammps_to_xyz_cli(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis lammps-to-xyz",
        description="Convert LAMMPS dump to multi-frame XYZ.",
    )
    parser.add_argument("input", help="LAMMPS .lammpstrj file")
    parser.add_argument("output", help="Output .xyz file")
    args = parser.parse_args(argv)
    lammps_to_xyz(args.input, args.output)


def _extract_frames_cli(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis extract-frames",
        description="Extract specific frames from an XYZ trajectory.",
    )
    parser.add_argument("input", help="Input XYZ trajectory")
    parser.add_argument("output", nargs="?", default="extracted.xyz",
                        help="Output XYZ file (default: extracted.xyz)")
    parser.add_argument(
        "--frames", default=None,
        help=(
            "Frame spec, e.g. '0-5,10:50:5,100:500:50'. "
            "Default: built-in logarithmic set."
        ),
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-c", "--count", action="store_true",
                        help="Only count frames, do not extract")
    args = parser.parse_args(argv)

    if args.count:
        n = count_xyz_frames_quick(args.input)
        print(f"Total frames: {n}")
        return

    keep = parse_frame_spec(args.frames) if args.frames else None
    extract_frames(args.input, args.output, keep, verbose=args.verbose)
