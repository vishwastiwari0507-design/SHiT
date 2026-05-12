"""
SHiT
====
A Python package for molecular dynamics trajectory analysis.

Modules
-------
core              – shared trajectory I/O, geometry, parsers
analysis          – bond order, RDF, dissociation, water density,
                    surface coverage, atom density
convert           – format conversion (LAMMPS→XYZ) and frame extraction

Quick start
-----------
>>> from SHiT.analysis import run_rdf
>>> results = run_rdf("traj.xyz", "input.txt", output_dir="rdf_out")

CLI
---
    shit --help
"""

from .analysis import (
    run_bond_order,
    run_rdf,
    run_dissociation,
    run_water_density,
    run_surface_coverage,
    run_atom_density_z,
)
from .convert import lammps_to_xyz, extract_frames

__version__ = "1.0.0"
__all__ = [
    "run_bond_order",
    "run_rdf",
    "run_dissociation",
    "run_water_density",
    "run_surface_coverage",
    "run_atom_density_z",
    "lammps_to_xyz",
    "extract_frames",
]
