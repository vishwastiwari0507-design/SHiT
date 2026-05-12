"""
mdanalysis.core.trajectory
==========================
Shared trajectory reading utilities for XYZ and LAMMPS formats.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Generator, List, Set, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class Frame:
    """One snapshot from a trajectory."""
    index: int
    elements: List[str]
    coords: np.ndarray          # shape (n_atoms, 3)

    def __post_init__(self):
        self.coords = np.asarray(self.coords, dtype=float)

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    def filter_elements(self, keep: Set[str]) -> "Frame":
        """Return a new Frame keeping only atoms whose element is in *keep*."""
        mask = [e in keep for e in self.elements]
        return Frame(
            index=self.index,
            elements=[e for e, m in zip(self.elements, mask) if m],
            coords=self.coords[mask],
        )

    def atoms_of(self, element: str) -> Tuple[List[int], np.ndarray]:
        """Return (local_indices, positions) for all atoms of *element*."""
        idx = [i for i, e in enumerate(self.elements) if e == element]
        pos = self.coords[idx] if idx else np.empty((0, 3))
        return idx, pos


@dataclass
class LammpsFrame:
    """One snapshot from a LAMMPS trajectory (lammpstrj)."""
    timestep: int
    n_atoms: int
    box_bounds: List[List[float]]   # [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
    atom_types: List[int]
    atom_data: dict                 # column_name -> list of values


# ---------------------------------------------------------------------------
# XYZ reader
# ---------------------------------------------------------------------------

def read_xyz(path: str) -> Generator[Frame, None, None]:
    """Yield Frame objects from a (multi-frame) XYZ file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"XYZ file not found: {path!r}")

    with open(path) as fh:
        frame_idx = 0
        while True:
            line = fh.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            try:
                n_atoms = int(line)
            except ValueError:
                raise ValueError(f"Expected atom count, got: {line!r}")

            fh.readline()           # comment line – skip

            elements, coords = [], []
            for _ in range(n_atoms):
                parts = fh.readline().split()
                if len(parts) < 4:
                    raise ValueError("Atom line must have: element x y z")
                elements.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

            yield Frame(index=frame_idx, elements=elements, coords=np.array(coords))
            frame_idx += 1


def load_xyz(path: str) -> List[Frame]:
    """Load all frames from an XYZ file into a list."""
    return list(read_xyz(path))


def count_xyz_frames(path: str) -> int:
    """Count frames in an XYZ file without loading coordinates."""
    n = 0
    with open(path) as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            try:
                n_atoms = int(line.strip())
            except ValueError:
                continue
            fh.readline()
            for _ in range(n_atoms):
                fh.readline()
            n += 1
    return n


# ---------------------------------------------------------------------------
# LAMMPS reader
# ---------------------------------------------------------------------------

def read_lammpstrj(path: str) -> Generator[LammpsFrame, None, None]:
    """Yield LammpsFrame objects from a LAMMPS dump file."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"LAMMPS trajectory not found: {path!r}")

    with open(path) as fh:
        lines = fh.readlines()

    i = 0
    total = len(lines)

    while i < total:
        if not lines[i].startswith("ITEM: TIMESTEP"):
            i += 1
            continue

        i += 1
        timestep = int(lines[i].strip()); i += 1

        while i < total and "ITEM: NUMBER OF ATOMS" not in lines[i]:
            i += 1
        i += 1
        n_atoms = int(lines[i].strip()); i += 1

        while i < total and "ITEM: BOX BOUNDS" not in lines[i]:
            i += 1
        i += 1
        box_bounds = []
        for _ in range(3):
            box_bounds.append([float(v) for v in lines[i].strip().split()])
            i += 1

        while i < total and "ITEM: ATOMS" not in lines[i]:
            i += 1
        header_parts = lines[i].strip().split()[2:]   # drop "ITEM: ATOMS"
        col_idx = {col: j for j, col in enumerate(header_parts)}
        i += 1

        atom_types: List[int] = []
        raw_data: dict = {col: [] for col in header_parts}

        for _ in range(n_atoms):
            parts = lines[i].strip().split(); i += 1
            for col, j in col_idx.items():
                raw_data[col].append(parts[j])
            atom_types.append(int(parts[col_idx["type"]]))

        yield LammpsFrame(
            timestep=timestep,
            n_atoms=n_atoms,
            box_bounds=box_bounds,
            atom_types=atom_types,
            atom_data=raw_data,
        )
