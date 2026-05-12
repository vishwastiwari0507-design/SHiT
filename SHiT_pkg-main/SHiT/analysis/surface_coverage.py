"""
mdanalysis.analysis.surface_coverage
======================================
Surface coverage analysis along the z-direction.

CLI
---
    mdanalysis surface-coverage -d traj.xyz -i input.txt -p para.txt [-o output.txt]

input.txt::

    Adsorbent:
    Si

    Adsorbates:
    O
    H
    OH

para.txt::

    Surface_Parameters:
    Adsorbent_buffer = 1.5
    Positive_z_range = 5.0
    Negative_z_range = 1.0

    Adsorbent-Adsorbate:

    Si-H
    Bondlength:
    Si-H = 1.1-1.8

    Si-O
    Bondlength:
    Si-O = 1.4-1.7

    Si-OH
    Bondlength:
    Si-O = 1.3-1.7
    O-H = 0.8-1.2
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

import numpy as np

from ..core.trajectory import load_xyz
from ..core.parsers import parse_range


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SurfaceParams:
    buffer: float = 2.0
    pos_z: float = 3.0
    neg_z: float = 2.0


@dataclass
class BondCriteria:
    e1: str
    e2: str
    lo: float
    hi: float

    def ok(self, d: float) -> bool:
        return self.lo <= d <= self.hi


@dataclass
class AdsDef:
    name: str
    criteria: List[BondCriteria]

    @property
    def is_bridge(self):
        return len(self.criteria) > 1


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

def _parse_input(path: str) -> Tuple[List[str], List[str]]:
    adsorbents, adsorbates = [], []
    section = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            low = line.lower()
            if low.startswith("adsorbent:"):
                section = "ads"; continue
            if low.startswith("adsorbates:"):
                section = "adsorbates"; continue
            if section == "ads":
                adsorbents.append(line)
            elif section == "adsorbates":
                adsorbates.append(line)
    return adsorbents, adsorbates


def _parse_para(path: str) -> Tuple[SurfaceParams, Dict[str, AdsDef]]:
    sp = SurfaceParams()
    defs: Dict[str, AdsDef] = {}

    with open(path) as fh:
        content = fh.read()

    # Surface parameters
    for attr, pattern in [
        ("buffer", r"Adsorbent_buffer\s*=\s*([\d.]+)"),
        ("pos_z",  r"Positive_z_range\s*=\s*([\d.]+)"),
        ("neg_z",  r"Negative_z_range\s*=\s*([\d.]+)"),
    ]:
        m = re.search(pattern, content, re.IGNORECASE)
        if m:
            setattr(sp, attr, float(m.group(1)))

    # Adsorption definitions (after "Adsorbent-Adsorbate:")
    parts = re.split(r"Adsorbent-Adsorbate:", content, flags=re.IGNORECASE)
    if len(parts) < 2:
        return sp, defs

    ads_section = parts[1]
    current_name = None
    current_criteria: List[BondCriteria] = []

    def _flush():
        nonlocal current_name, current_criteria
        if current_name and current_criteria:
            defs[current_name] = AdsDef(name=current_name, criteria=list(current_criteria))
        current_name = None
        current_criteria = []

    for line in ads_section.split("\n"):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        low = line.lower()
        if low.startswith("bondlength") or low.startswith("bond_length") \
                or low.startswith("bond length"):
            continue
        # Bond spec: X-Y = lo-hi
        bm = re.match(r"([A-Za-z]+)-([A-Za-z]+)\s*=\s*([\d.]+)-([\d.]+)", line)
        if bm:
            current_criteria.append(BondCriteria(
                e1=bm.group(1), e2=bm.group(2),
                lo=float(bm.group(3)), hi=float(bm.group(4)),
            ))
            continue
        # Pattern name: X-Y (no digits, no =)
        pm = re.match(r"^([A-Za-z]+-[A-Za-z]+)$", line)
        if pm:
            _flush()
            current_name = pm.group(1)
            continue

    _flush()
    return sp, defs


# ---------------------------------------------------------------------------
# Surface identification
# ---------------------------------------------------------------------------

def _identify_surface(atoms_z: np.ndarray, adsorbent_mask: np.ndarray,
                      sp: SurfaceParams) -> Tuple[float, float, float]:
    """Return (z_avg, z_lower, z_upper) for the active surface."""
    z_ads = atoms_z[adsorbent_mask]
    z_max = z_ads.max()
    near = z_ads[z_ads >= z_max - sp.buffer]
    z_avg = near.mean()
    return z_avg, z_avg - sp.neg_z, z_avg + sp.pos_z


# ---------------------------------------------------------------------------
# Adsorption detection
# ---------------------------------------------------------------------------

def _euclidean(p1, p2) -> float:
    return float(np.linalg.norm(p1 - p2))


def _detect_simple(ads_def: AdsDef, surf_mask: np.ndarray,
                   region_mask: np.ndarray,
                   elements: List[str], coords: np.ndarray,
                   excluded: Set[int]) -> List[int]:
    """Return list of adsorbate atom indices that are bonded to any surface atom."""
    crit = ads_def.criteria[0]
    result = []
    region_idx = [i for i in range(len(elements))
                  if region_mask[i] and elements[i] == crit.e2
                  and i not in excluded]
    surf_idx = [i for i in range(len(elements))
                if surf_mask[i] and elements[i] == crit.e1]
    for ri in region_idx:
        for si in surf_idx:
            if crit.ok(_euclidean(coords[ri], coords[si])):
                result.append(ri)
                break
    return result


def _detect_bridge(ads_def: AdsDef, surf_mask: np.ndarray,
                   region_mask: np.ndarray,
                   elements: List[str], coords: np.ndarray,
                   excluded: Set[int]) -> List[int]:
    """Return list of terminal adsorbate indices in A-C-B bridge."""
    c1, c2 = ads_def.criteria[0], ads_def.criteria[1]
    result = []
    mid_idx = [i for i in range(len(elements))
               if region_mask[i] and elements[i] == c1.e2 and i not in excluded]
    term_idx = [i for i in range(len(elements))
                if region_mask[i] and elements[i] == c2.e2 and i not in excluded]
    surf_idx = [i for i in range(len(elements))
                if surf_mask[i] and elements[i] == c1.e1]

    used_mid: Set[int] = set()
    used_term: Set[int] = set()

    for mi in mid_idx:
        if mi in used_mid:
            continue
        bonded_surf = any(c1.ok(_euclidean(coords[mi], coords[si])) for si in surf_idx)
        if not bonded_surf:
            continue
        for ti in term_idx:
            if ti in used_term or ti == mi:
                continue
            if c2.ok(_euclidean(coords[mi], coords[ti])):
                used_mid.add(mi); used_term.add(ti)
                result.append(mi)   # we count the mid-atom (e.g. O in Si-OH)
                break

    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def run_surface_coverage(
    trajectory: str,
    input_file: str,
    para_file: str,
    output_file: str = "surface_coverage.txt",
) -> List[Dict]:
    """Run surface coverage analysis.

    Returns list of per-frame dicts.
    """
    adsorbents, adsorbates = _parse_input(input_file)
    sp, ads_defs = _parse_para(para_file)

    frames = load_xyz(trajectory)
    results = []
    n = len(frames)

    # Sort: bridge first
    sorted_defs = sorted(ads_defs.values(), key=lambda d: (not d.is_bridge, d.name))

    for frame in frames:
        elems = frame.elements
        coords = frame.coords
        z = coords[:, 2]

        ads_mask = np.array([e in adsorbents for e in elems])
        if not ads_mask.any():
            results.append(dict(frame=frame.index, counts={}, total=0))
            continue

        z_avg, z_lo, z_hi = _identify_surface(z, ads_mask, sp)

        surf_mask = ads_mask & (z >= z_avg - sp.buffer)
        region_mask = (z >= z_lo) & (z <= z_hi)

        excluded: Set[int] = set()
        counts: Dict[str, int] = {}

        for ads_def in sorted_defs:
            if ads_def.is_bridge:
                matched = _detect_bridge(
                    ads_def, surf_mask, region_mask, elems, coords, excluded
                )
            else:
                matched = _detect_simple(
                    ads_def, surf_mask, region_mask, elems, coords, excluded
                )
            counts[ads_def.name] = len(matched)
            excluded.update(matched)

        total = sum(counts.values())
        results.append(dict(frame=frame.index, counts=counts, total=total))

        if frame.index % 100 == 0 or frame.index == n - 1:
            print(f"  Frame {frame.index+1}/{n}: total={total}  {counts}")

    # Write output
    all_types = sorted({k for r in results for k in r["counts"]})
    with open(output_file, "w") as fh:
        fh.write("Frame\t" + "\t".join(all_types) + "\tTotal\n")
        fh.write("-" * (8 + 8 * len(all_types) + 8) + "\n")
        for r in results:
            row = f"{r['frame']}"
            for t in all_types:
                row += f"\t{r['counts'].get(t, 0)}"
            row += f"\t{r['total']}\n"
            fh.write(row)
        fh.write("\n# Summary\n")
        for t in all_types:
            vals = [r["counts"].get(t, 0) for r in results]
            fh.write(f"# {t}: mean={np.mean(vals):.3f} std={np.std(vals):.3f}\n")
        totals = [r["total"] for r in results]
        fh.write(f"# Total: mean={np.mean(totals):.3f} std={np.std(totals):.3f}\n")

    print(f"\nResults written to {output_file}")
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="mdanalysis surface-coverage",
        description="Surface coverage analysis along z-direction.",
    )
    parser.add_argument("-d", "--trajectory", required=True)
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-p", "--parameters", required=True)
    parser.add_argument("-o", "--output", default="surface_coverage.txt")
    args = parser.parse_args(argv)

    run_surface_coverage(
        args.trajectory, args.input, args.parameters, args.output
    )


if __name__ == "__main__":
    main()
