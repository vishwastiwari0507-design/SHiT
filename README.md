# <img src="https://raw.githubusercontent.com/vishwastiwari0507-design/SHiT/main/SHiT_pkg-main/logo%20copy.png" alt="SHiT" width="25%">


![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/vishwastiwari0507-design/SHiT)
[![Paper](https://img.shields.io/badge/Paper-Publication-blue)](YOUR_PAPER_LINK)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)




## Table of contents
- [Description](#description)
- [Supported calculation types and simulation engines](#supported-calculations-and-simulation-engines)
- [Installation](#installation)
  - [(i) Installation from source](#i-installation-from-source)
  - [(ii) Installation via pip](#ii-installation-via-pip)
- [Usage](#usage)
  - [Command line interface](#command-line-interface)
  - [Python API](#python-api)
- [Examples](#examples)
- [Citation](#citation)
- [License](#license)



## Installation


### (i) Installation from source 

1. Create a virtual environment:
   ```bash
   conda create -n sea 
   conda activate sea
   ```

2. Install SEA_toolkit:
   ```bash
   git clone https://github.com/vishwastiwari0507-design/SEA_toolkit.git
   cd SEA_toolkit.git
   pip install -r requirements.txt
   pip install .
   ```

### (ii) Installation via pip

```bash
conda create -n sea python=3.12
conda activate sea
pip install SEA_toolkit
```

## Usage

### Command line interface

All tools are available under a single `shit` command:

```
sea <command> [options]
```

### Available commands

| Command | Original script | Description |
|---|---|---|
| `bond-order` | `bond_analysis.py` | Count di/tri/four-atom patterns |
| `rdf` | `rdf.py` | Radial distribution functions |
| `dissociation` | `Dissociation_H2O.py` | Water species classification |
| `water-density` | `Water_density_analysis.py` | Z-directional water density |
| `surface-coverage` | `Surface.Coverage.py` | Surface adsorbate coverage |
| `atom-density` | `atom_vs_z_direction.py` | Atom density along Z from XYZ |
| `lammps-to-xyz` | `lammpstrj_to_xyz.py` | Convert LAMMPS dump → XYZ |
| `extract-frames` | `pp.py` / `extracting_specific_frames.py` | Extract frames from XYZ |

---

### Python API

All tools are also importable as functions:

```python
from SHiT.analysis import (
    run_bond_order,
    run_rdf,
    run_dissociation,
    run_water_density,
    run_surface_coverage,
    run_atom_density_z,
)
from SHiT.convert import lammps_to_xyz, extract_frames

# Bond order
results = run_bond_order("traj.xyz", "input.txt", "para.txt", output_dir="out/")

# RDF
rdfs = run_rdf("traj.xyz", "input.txt", output_dir="rdf/")
r, g_r = rdfs["O-O"]   # numpy arrays

# Dissociation
frames = run_dissociation("traj.xyz", "input.txt", "para.txt", output_file="out.txt")

# Water density
bin_centers, density_matrix = run_water_density("traj.xyz", "input.txt", "para.txt")

# Surface coverage
cov = run_surface_coverage("traj.xyz", "input.txt", "para.txt")

# Atom density Z
hists = run_atom_density_z("traj.xyz", "input.txt", "para.txt")

# Conversion
lammps_to_xyz("dump.lammpstrj", "out.xyz")
extract_frames("traj.xyz", "subset.xyz", frames_to_keep={0,1,2,100,500})
```

---
