"""
SHiT – unified command-line dispatcher.

Usage
-----
    shit <command> [options]

Commands
--------
    bond-order        Count molecular patterns (di/tri/four-atom)
    rdf               Radial distribution functions
    dissociation      Water dissociation species (H2O, OH-, H3O+, …)
    water-density     Z-directional water density profile
    surface-coverage  Adsorbate surface coverage along z
    atom-density      Atom density along Z from XYZ trajectory
    lammps-to-xyz     Convert LAMMPS dump → multi-frame XYZ
    extract-frames    Extract specific frames from XYZ trajectory
"""

import sys


COMMANDS = {
    "bond-order":       "SHiT.analysis.bond_order:main",
    "rdf":              "SHiT.analysis.rdf:main",
    "dissociation":     "SHiT.analysis.dissociation:main",
    "water-density":    "SHiT.analysis.water_density:main",
    "surface-coverage": "SHiT.analysis.surface_coverage:main",
    "atom-density":     "SHiT.analysis.atom_density_z:main",
    "lammps-to-xyz":    "SHiT.convert:_lammps_to_xyz_cli",
    "extract-frames":   "SHiT.convert:_extract_frames_cli",
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print(__doc__)
        print("Available commands:")
        for cmd in COMMANDS:
            print(f"    {cmd}")
        sys.exit(0)

    cmd = sys.argv[1]
    if cmd not in COMMANDS:
        print(f"Unknown command: {cmd!r}\n")
        print("Run 'shit --help' for available commands.")
        sys.exit(1)

    module_path, func_name = COMMANDS[cmd].rsplit(":", 1)
    import importlib
    module = importlib.import_module(module_path)
    func = getattr(module, func_name)
    func(sys.argv[2:])


if __name__ == "__main__":
    main()
