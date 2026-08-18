
# <img src=https://github.com/vishwastiwari0507-design/SHiT/blob/main/SHiT_pkg-main/logo2.png alt="SHiT" width="25%">

![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/vishwastiwari0507-design/SHiT)
[![Paper](https://img.shields.io/badge/Paper-Publication-blue)](YOUR_PAPER_LINK)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)




## Table of contents
- [Description](#description)
- [Installation](#installation)
  - [(i) Installation from source](#i-installation-from-source)
  - [(ii) Installation via pip](#ii-installation-via-pip)
- [Usage](#usage)
  - [Command line interface](#command-line-interface)
  - [Python API](#python-api)
- [Examples](#Examples)
- [Package structure](#Package-structure)
- [Citation](#citation)

## Description
SEA-toolkit (Solid Electrolyte Analysis) is a Python package for post-processing molecular dynamics trajectories. It provides a unified CLI and API for bond-order counting, RDF, water dissociation, water/atom density profiles, and surface coverage analysis of solid electrolyte interfaces.

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
from sea.analysis import (
    run_bond_order,
    run_rdf,
    run_dissociation,
    run_water_density,
    run_surface_coverage,
    run_atom_density_z,
)
from sea.convert import lammps_to_xyz, extract_frames

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


## Examples

### Bond-order analysis
```bash
sea bond-order -d traj.xyz -i input.txt -p para.txt
sea bond-order -d traj.xyz -i input.txt -p para.txt -o results/
```

**input.txt** – one pattern per line:
```
Si-O
Si-O-Si
Si-H
```

**para.txt** – bond / angle / dihedral ranges:
```
Si-O
 bond_length:
   Si-O = 1.5-1.6

Si-H
 bond_length:
   Si-H = 1.5-1.8

Si-O-Si
 bond_length:
   Si-O = 1.5-1.7
   O-Si = 1.5-1.7
 bond_angle:
   Si-O-Si = 130-160
```

Output: one `<pattern>.txt` per pattern with `frame\tcount`.

---

### RDF analysis
```bash
sea rdf -d traj.xyz -i input.txt
sea rdf -d traj.xyz -i input.txt -o rdf_out/
```

**input.txt**:
```
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
```

Output: `O-O.txt`, `O-H.txt`, `Si-O.txt`, and `all.txt`.

---

### Water dissociation
```bash
sea dissociation -d traj.xyz -i input.txt -p para.txt
sea dissociation -d traj.xyz -i input.txt -p para.txt -o results.txt -v
```

**input.txt**:
```
Dimension:
a = 0 38.241
b = 0 38.241
c = 0 57.32

Region:
A = 0 38.241
B = 0 38.241
C = 27 57.32

PBC:
T T F
```

**para.txt**:
```
Bond length:
H-O = 0.8 1.2
H-Si = 1.5 1.8
O-Si = 1.4 1.7

Bond angle:
H-O-H = 90 115
```

Output columns: `Frame  H+/H3O+  OH-  Total  H2O  Si-OH  Si-H`

---

### Water density profile
```bash
sea water-density -d traj.xyz -i input.txt -p para.txt
sea water-density -d traj.xyz -i input.txt -p para.txt -o density_out/
```

**input.txt**:
```
Dimension:
a = 38.2
b = 38.2
c = 57

PBC:
T T T

Bin_size = 2
```

**para.txt**:
```
Bond_length:
O-H = 0.8-1.2

Bond_angle:
H-O-H = 90-115
```

Outputs: `H2O.txt`, `density_profile.txt`, `water_count.txt`

---

### Surface coverage
```bash
sea surface-coverage -d traj.xyz -i input.txt -p para.txt
sea surface-coverage -d traj.xyz -i input.txt -p para.txt -o coverage.txt
```

**input.txt**:
```
Adsorbent:
Si

Adsorbates:
O
H
OH
```

**para.txt**:
```
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
```

---

### Atom density along Z
```bash
sea atom-density -d traj.xyz -i input.txt -p para.txt
sea atom-density -d traj.xyz -i input.txt -p para.txt -o z_plots/
```

**input.txt**:
```
Dimension:
a = 0:38.241
b = 0:38.241
c = 0:57.32

PBC:
T T F

Bin_size = 0.5
```

**para.txt**:
```
Elements:
Si
O
H
```

Outputs: `<element>_z_data.txt` and `<element>_z_distribution.png` per element, plus `ALL_elements_z_distribution.png`.

---

### Format conversion
```bash
# Convert LAMMPS dump to XYZ
sea lammps-to-xyz trajectory.lammpstrj trajectory.xyz

# Extract frames (default logarithmic set)
sea extract-frames trajectory.xyz extracted.xyz

# Extract custom frame set
sea extract-frames trajectory.xyz extracted.xyz --frames "0-5,10:50:5,100:500:50"

# Just count frames
sea extract-frames trajectory.xyz --count
```

---

### Input file format notes

#### Range values
All input / parameter files accept flexible range formats:
- `1.5-1.7`   dash-separated
- `1.5:1.7`   colon-separated
- `1.5 1.7`   space-separated

#### PBC flags
`T T F` = periodic in x and y, non-periodic in z.

---

## Package structure

```
SEA_pkg/
├── pyproject.toml
├── README.md
└── SEA/
    ├── __init__.py
    ├── cli.py                     ← unified CLI dispatcher  (shit <command>)
    ├── core/
    │   ├── trajectory.py          ← XYZ + LAMMPS readers
    │   ├── geometry.py            ← distances, angles, PBC
    │   └── parsers.py             ← input / para file parsers
    ├── analysis/
    │   ├── bond_order.py
    │   ├── rdf.py
    │   ├── dissociation.py
    │   ├── water_density.py
    │   ├── surface_coverage.py
    │   └── atom_density_z.py
    └── convert/
        └── __init__.py            ← lammps_to_xyz, extract_frames
```

## Citation

If you use SEA_toolkit in your research, please cite:

Tiwari, V. and Gu, C. and Elgabarty, H. and Brehm, M. "Atomistic Insights into Si(100) Surface Passivation by Water: ReaxFF Molecular Dynamics and Fine-Tuned Machine Learning Interatomic Potential Simulations." arXiv preprint arXiv:--- (2026).

```
@misc{------,
      title={}, 
      author={Vishwas Tiwari and Chenying Gu and Hossam Elgabarty and Martin Brehm},
      year={2026},
      eprint={-},
      archivePrefix={arXiv},
      primaryClass={-},
      url={https://arxiv.org/}, 
}
```
