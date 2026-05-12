"""
mdanalysis.core.parsers
=======================
Shared parsers for para.txt and input.txt files used across all analyses.
"""

from __future__ import annotations

import re
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# Range helpers
# ---------------------------------------------------------------------------

def parse_range(s: str) -> Tuple[float, float]:
    """Parse 'min-max', 'min:max', or 'min max' into (min, max).

    Handles negative numbers like '-5-20' or '-10--5'.
    """
    s = s.strip()

    # Try colon notation first: 'start:end'
    if ":" in s:
        parts = s.split(":", 1)
        return float(parts[0]), float(parts[1])

    # Try space notation: 'min max'
    parts = s.split()
    if len(parts) == 2:
        return float(parts[0]), float(parts[1])

    # Dash notation with possible negatives
    m = re.match(r"^([+-]?\d+(?:\.\d*)?)[ \t]*-[ \t]*([+-]?\d+(?:\.\d*)?)$", s)
    if m:
        return float(m.group(1)), float(m.group(2))

    # Fall back: extract all numbers
    nums = re.findall(r"-?\d+\.?\d*", s)
    if len(nums) >= 2:
        lo, hi = float(nums[0]), float(nums[1])
        return (lo, hi) if lo <= hi else (hi, lo)
    if len(nums) == 1:
        v = float(nums[0])
        return (0.0, v) if v >= 0 else (v, 0.0)

    raise ValueError(f"Cannot parse range: {s!r}")


# ---------------------------------------------------------------------------
# Bond-order / surface-coverage parameter file
# ---------------------------------------------------------------------------

def parse_bond_para(path: str) -> Dict:
    """Parse a bond-parameter file used by bond_order and surface_coverage.

    Returns
    -------
    dict keyed by pattern string e.g. 'Si-O', 'Si-O-Si', …
    Each value is a dict with keys 'bond_length', 'bond_angle', 'dihedral'.
    """
    params: Dict = {}
    current_pattern = None
    current_section = None

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            if ":" not in line and "=" not in line:
                current_pattern = line
                params.setdefault(current_pattern, {})
                current_section = None
                continue

            if line.endswith(":") and "=" not in line:
                key = line[:-1].strip().lower().replace(" ", "_")
                if key in ("bond_length", "bondlength", "bond_angle",
                           "bondangle", "dihedral"):
                    # Normalise
                    if key in ("bondlength",):
                        key = "bond_length"
                    elif key in ("bondangle",):
                        key = "bond_angle"
                    current_section = key
                    if current_pattern:
                        params[current_pattern].setdefault(current_section, {})
                else:
                    current_section = None
                continue

            if "=" in line and current_pattern and current_section:
                lhs, rhs = line.split("=", 1)
                label = lhs.strip()
                try:
                    lo, hi = parse_range(rhs.strip())
                    params[current_pattern][current_section][label] = (lo, hi)
                except ValueError:
                    pass

    return params


# ---------------------------------------------------------------------------
# Box / region / PBC input file (Dissociation, RDF, water density)
# ---------------------------------------------------------------------------

class BoxInputParser:
    """Parse a simulation-box input file.

    Supports several value formats:
    - 'a = 0-20'      (dash separated)
    - 'a = 0:20'      (colon separated)
    - 'a = 0 20'      (space separated)
    - 'a = 20'        (single value → 0..value)
    """

    def __init__(self, path: str):
        self.path = path
        self.box_min = [0.0, 0.0, 0.0]
        self.box_max = [0.0, 0.0, 0.0]
        self.region_min = [0.0, 0.0, 0.0]
        self.region_max = [0.0, 0.0, 0.0]
        self.pbc = [True, True, False]
        self._region_set = False
        self._extra: Dict = {}

    # ------------------------------------------------------------------
    def parse(self) -> "BoxInputParser":
        section = None
        with open(self.path) as fh:
            for lineno, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue

                low = line.lower()
                if low.startswith("dimension"):
                    section = "dim"; continue
                if low.startswith("region"):
                    section = "region"; self._region_set = True; continue
                if low.startswith("pbc"):
                    section = "pbc"; continue
                if "=" in line and not any(
                    line.lower().startswith(k)
                    for k in ("dimension", "region", "pbc")
                ):
                    self._parse_kv(line, section, lineno)
                    continue
                if section == "pbc":
                    parts = line.split()
                    if len(parts) >= 3:
                        self.pbc = [p.upper() == "T" for p in parts[:3]]

        if not self._region_set:
            self.region_min = list(self.box_min)
            self.region_max = list(self.box_max)

        return self

    # ------------------------------------------------------------------
    _DIM_MAP = {"a": 0, "b": 1, "c": 2}
    _REG_MAP = {"a": 0, "b": 1, "c": 2}

    # Keys that are always plain scalars, never ranges
    _SCALAR_KEYS = {
        "r_max", "rmax", "r_cutoff", "cutoff",
        "bin_size", "binsize", "bin_width", "binwidth", "dr",
    }

    def _parse_kv(self, line: str, section, lineno: int):
        if "=" not in line:
            return
        key, val = line.split("=", 1)
        key = key.strip().lower()
        val = val.strip()

        # Always try plain float first for known scalar keys
        if key in self._SCALAR_KEYS:
            try:
                self._extra[key] = float(val)
                return
            except ValueError:
                pass

        # Try to parse as a range
        try:
            lo, hi = parse_range(val)
        except ValueError:
            # Fall back: try plain float for unknown keys
            try:
                self._extra[key] = float(val)
            except ValueError:
                self._extra[key] = val
            return

        if section == "dim" and key in self._DIM_MAP:
            i = self._DIM_MAP[key]
            self.box_min[i] = lo
            self.box_max[i] = hi
        elif section == "region" and key in self._REG_MAP:
            i = self._REG_MAP[key]
            self.region_min[i] = lo
            self.region_max[i] = hi
        else:
            self._extra[key] = (lo, hi)

    # ------------------------------------------------------------------
    @property
    def box_dims(self):
        import numpy as np
        return np.array(self.box_max) - np.array(self.box_min)

    @property
    def region_dims(self):
        import numpy as np
        return np.array(self.region_max) - np.array(self.region_min)

    def extra(self, key: str, default=None):
        return self._extra.get(key, default)
