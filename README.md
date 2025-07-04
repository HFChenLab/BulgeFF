# BulgeFF: Bulge-Specific Force Field for Accurate RNA Sampling

## Overview

BulgeFF is a bulge-specific force field that applies real-time energy corrections to the η-θ pseudo-torsions of target bulge nucleotides during MD simulations. By guiding proper conformational sampling, BulgeFF significantly improves the simulation accuracy of RNA bulges. Users can use BulgeFF by simply provide the initial structure and essential input information to rapidly generate PLUMED input files for MD simulations.

## Installation

```bash
git clone https://github.com/yourusername/BulgeFF.git
cd BulgeFF
pip install -r requirements.txt
```

## Requirements

- Python 3.7+
- NumPy
- BioPython
- MDAnalysis
- PLUMED (for MD simulation)
- GROMACS (for MD simulations)


## Quick Start

### Basic Usage

```bash
python BulgeFF.py \
    --bulge_pdb 2jym.pdb \
    --md_pdb reference.pdb \
    --bulge_name G \
    --bulge_id 6
```

### Command Line Arguments

```bash
python BulgeFF.py [OPTIONS]

Options:
  --bulge_pdb TEXT     Name of PDB file containing bulge structure [default: 2jym.pdb]
  --md_pdb TEXT        Name of PDB file for MD simulation [default: reference.pdb]
  --bulge_name TEXT... List of bulge residue names (A, U, G, C) [default: G]
  --bulge_id INT...    List of bulge residue IDs [default: 6]
  --help              Show this message and exit
```

### Usage Examples

**Single bulge residue:**
```bash
python BulgeFF.py --bulge_pdb 2jym.pdb --md_pdb reference.pdb --bulge_name G --bulge_id 6
```

**Multiple bulge residues:**
```bash
python BulgeFF.py \
    --bulge_pdb rna_structure.pdb \
    --md_pdb md_reference.pdb \
    --bulge_name G A U C \
    --bulge_id 25 26 36 37
```

## Input Requirements

### Required Files
- **Bulge PDB File**: Input RNA structure containing bulge regions (e.g., `2jym.pdb`)
- **MD PDB File**: Structure file obtained after energy minimization with GROMACS (ensures atom indexing consistency) (e.g., `reference.pdb`)
- **Bulge Specification**: 
  - Single-letter names of bulge nucleotides (A, U, G, C)
  - Corresponding residue IDs in the structure

## Output Files

The program generates `plumed.dat` - the PLUMED input file with energy correction terms for GROMACS simulations.

## Integration with GROMACS

The generated PLUMED files can be directly used with GROMACS:

```bash
gmx mdrun -deffnm md -plumed plumed.dat
```
