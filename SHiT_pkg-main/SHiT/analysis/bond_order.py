"""
mdanalysis.analysis.bond_order
==============================
Count diatomic, triatomic, and four-atom molecular patterns in XYZ trajectories.

CLI
---
    mdanalysis bond-order -d traj.xyz -i input.txt -p para.txt

input.txt  – one pattern per line, e.g.::

    Si-O
    Si-O-Si
    Si-H

para.txt   – bond / angle / dihedral ranges, e.g.::

    Si-O
     bond_length:
       Si-O = 1.5-1.6

Output
------
One text file per pattern: ``<pattern>.txt``  with columns ``frame\\tcount``.
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Tuple

from ..core.trajectory import Frame, load_xyz
from ..core.parsers import parse_bond_para
from ..core.geometry import distance_3d, angle_3d, dihedral_3d


# ---------------------------------------------------------------------------
# Detection helpers
# ---------------------------------------------------------------------------

def _bond_range(params: Dict, pattern: str, a: str, b: str) -> Tuple[float, float]:
    bl = params[pattern].get("bond_length", {})
    for key in (f"{a}-{b}", f"{b}-{a}"):
        if key in bl:
            return bl[key]
    raise ValueError(
        f"No bond_length for {a}-{b} in pattern '{pattern}' inside para file."
    )


def detect_diatomic(pattern: str, frames: List[Frame], params: Dict) -> List[int]:
    e1, e2 = pattern.split("-")
    lo, hi = _bond_range(params, pattern, e1, e2)
    counts = []
    for frame in frames:
        idx1 = [i for i, e in enumerate(frame.elements) if e == e1]
        idx2 = [i for i, e in enumerate(frame.elements) if e == e2]
        found: set = set()
        same = e1 == e2
        for i in idx1:
            for j in idx2:
                if same and j <= i:
                    continue
                d = distance_3d(frame.coords[i], frame.coords[j])
                if lo <= d <= hi:
                    found.add(tuple(sorted((i, j))))
        counts.append(len(found))
    return counts


def detect_triatomic(pattern: str, frames: List[Frame], params: Dict) -> List[int]:
    eA, eB, eC = pattern.split("-")
    rAB = _bond_range(params, pattern, eA, eB)
    rBC = _bond_range(params, pattern, eB, eC)
    ba = params[pattern].get("bond_angle", {})
    ak = f"{eA}-{eB}-{eC}"
    if ak not in ba:
        raise ValueError(f"No bond_angle '{ak}' in pattern '{pattern}'.")
    ang_lo, ang_hi = ba[ak]

    counts = []
    for frame in frames:
        idxA = [i for i, e in enumerate(frame.elements) if e == eA]
        idxB = [i for i, e in enumerate(frame.elements) if e == eB]
        idxC = [i for i, e in enumerate(frame.elements) if e == eC]
        found: set = set()
        for iB in idxB:
            for iA in idxA:
                if iA == iB:
                    continue
                if not (rAB[0] <= distance_3d(frame.coords[iA], frame.coords[iB]) <= rAB[1]):
                    continue
                for iC in idxC:
                    if iC in (iA, iB):
                        continue
                    if not (rBC[0] <= distance_3d(frame.coords[iB], frame.coords[iC]) <= rBC[1]):
                        continue
                    ang = angle_3d(frame.coords[iA], frame.coords[iB], frame.coords[iC])
                    if ang_lo <= ang <= ang_hi:
                        found.add(tuple(sorted((iA, iB, iC))))
        counts.append(len(found))
    return counts


def detect_four_atom(pattern: str, frames: List[Frame], params: Dict) -> List[int]:
    eA, eB, eC, eD = pattern.split("-")
    rAB = _bond_range(params, pattern, eA, eB)
    rBC = _bond_range(params, pattern, eB, eC)
    rCD = _bond_range(params, pattern, eC, eD)
    ba = params[pattern].get("bond_angle", {})
    di = params[pattern].get("dihedral", {})
    ak1, ak2 = f"{eA}-{eB}-{eC}", f"{eB}-{eC}-{eD}"
    dk = f"{eA}-{eB}-{eC}-{eD}"
    for k in (ak1, ak2, dk):
        if k not in {**ba, **di}:
            raise ValueError(f"Missing parameter '{k}' for pattern '{pattern}'.")
    a1_lo, a1_hi = ba[ak1]
    a2_lo, a2_hi = ba[ak2]
    d_lo, d_hi = di[dk]

    counts = []
    for frame in frames:
        idxA = [i for i, e in enumerate(frame.elements) if e == eA]
        idxB = [i for i, e in enumerate(frame.elements) if e == eB]
        idxC = [i for i, e in enumerate(frame.elements) if e == eC]
        idxD = [i for i, e in enumerate(frame.elements) if e == eD]
        found: set = set()
        for iB in idxB:
            for iC in idxC:
                if iC == iB:
                    continue
                if not (rBC[0] <= distance_3d(frame.coords[iB], frame.coords[iC]) <= rBC[1]):
                    continue
                for iA in idxA:
                    if iA in (iB, iC):
                        continue
                    if not (rAB[0] <= distance_3d(frame.coords[iA], frame.coords[iB]) <= rAB[1]):
                        continue
                    ang1 = angle_3d(frame.coords[iA], frame.coords[iB], frame.coords[iC])
                    if not (a1_lo <= ang1 <= a1_hi):
                        continue
                    for iD in idxD:
                        if iD in (iA, iB, iC):
                            continue
                        if not (rCD[0] <= distance_3d(frame.coords[iC], frame.coords[iD]) <= rCD[1]):
                            continue
                        ang2 = angle_3d(frame.coords[iB], frame.coords[iC], frame.coords[iD])
                        if not (a2_lo <= ang2 <= a2_hi):
                            continue
                        dih = dihedral_3d(
                            frame.coords[iA], frame.coords[iB],
                            frame.coords[iC], frame.coords[iD],
                        )
                        if d_lo <= dih <= d_hi:
                            found.add(tuple(sorted((iA, iB, iC, iD))))
        counts.append(len(found))
    return counts


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_bond_order(
    trajectory: str,
    input_file: str,
    para_file: str,
    output_dir: str = ".",
) -> Dict[str, List[int]]:
    """Run bond-order analysis.

    Parameters
    ----------
    trajectory  : path to multi-frame XYZ file
    input_file  : path to input.txt listing patterns
    para_file   : path to para.txt with bond/angle ranges
    output_dir  : directory for output text files

    Returns
    -------
    dict mapping each pattern string to its per-frame count list.
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    patterns = _read_patterns(input_file)
    params = parse_bond_para(para_file)
    frames = load_xyz(trajectory)

    required = {e for pat in patterns for e in pat.split("-")}
    frames = [f.filter_elements(required) for f in frames]

    results: Dict[str, List[int]] = {}
    for pat in patterns:
        n = len(pat.split("-"))
        if pat not in params:
            raise ValueError(f"Pattern '{pat}' not found in para file.")
        if n == 2:
            counts = detect_diatomic(pat, frames, params)
        elif n == 3:
            counts = detect_triatomic(pat, frames, params)
        elif n == 4:
            counts = detect_four_atom(pat, frames, params)
        else:
            raise ValueError(f"Pattern '{pat}' has {n} atoms; only 2-4 supported.")

        out_path = os.path.join(output_dir, f"{pat}.txt")
        with open(out_path, "w") as fh:
            for frame, cnt in zip(frames, counts):
                fh.write(f"{frame.index}\t{cnt}\n")
        results[pat] = counts
        print(f"  {pat}: written to {out_path}")

    return results


def _read_patterns(path: str) -> List[str]:
    patterns = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                patterns.append(line)
    return patterns


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis bond-order",
        description="Count molecular patterns in an XYZ trajectory.",
    )
    parser.add_argument("-d", "--trajectory", required=True)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--parameters", required=True)
    parser.add_argument("-o", "--output-dir", default=".")
    args = parser.parse_args(argv)

    run_bond_order(args.trajectory, args.input, args.parameters, args.output_dir)


if __name__ == "__main__":
    main()
