#!/usr/bin/env python3
"""
Docking Preparation Script (no Meeko or MDAnalysis)
- Downloads PDB file
- Prepares protein structure with PDBFixer
- Converts receptor and ligands to PDBQT using Open Babel
- Ligand names are taken from the SDF title lines
- All output goes into docking_files/
"""

import os
import requests
from pdbfixer import PDBFixer
from openmm.app import PDBFile, Simulation, ForceField
from openmm import Platform, VerletIntegrator, unit
from rdkit import Chem
import pandas as pd

# Create necessary directories
os.makedirs("protein_structures", exist_ok=True)
os.makedirs("docking_files", exist_ok=True)

# 1. DOWNLOAD PDB STRUCTURE
pdb_id = "7LME"
pdb_path = f"protein_structures/{pdb_id}.pdb"
r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
with open(pdb_path, "w") as f:
    f.write(r.text)

# 2. PREPARE PROTEIN STRUCTURE
fixer = PDBFixer(filename=pdb_path)
forcefield = ForceField("amber/protein.ff14SB.xml")

fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(forcefield=forcefield)

system = forcefield.createSystem(fixer.topology)
integrator = VerletIntegrator(0.001 * unit.picoseconds)
platform = Platform.getPlatformByName("CPU")
simulation = Simulation(fixer.topology, system, integrator, platform)
simulation.context.setPositions(fixer.positions)

print("Minimizing energy...")
simulation.minimizeEnergy()
positions = simulation.context.getState(getPositions=True).getPositions()

minimized_pdb_path = f"docking_files/{pdb_id}_fixed.pdb"
with open(minimized_pdb_path, "w") as f:
    PDBFile.writeFile(fixer.topology, positions, f)

# 3. CONVERT PROTEIN TO PDBQT USING OPEN BABEL

# xh tells openbabel to keep the hydrogens
receptor_pdbqt_path = f"docking_files/{pdb_id}.pdbqt"
os.system(f"obabel -ipdb {minimized_pdb_path} -opdbqt -O {receptor_pdbqt_path} -xh")

# 4. PREPARE LIGANDS
df = pd.read_csv("US20240293380_examples.csv")
names = df["Name"].str.split().str[1:].str.join("_")
df["Compound Name"] = names
df[["SMILES", "Compound Name"]].to_csv("US20240293380_picked.smi", index=False, header=False, sep=" ")

# Generate 3D structures from SMILES with protonation
os.system("scrub.py US20240293380_picked.smi -o mols.sdf --ph_low 7.4 --ph_high 7.4")

# Ensure all Hs are present and write new SDF with clean names
supplier = Chem.SDMolSupplier("mols.sdf", removeHs=False)
writer = Chem.SDWriter("mols_with_H.sdf")

for mol in supplier:
    if mol is None:
        continue
    mol = Chem.AddHs(mol, addCoords=False)
    for prop in mol.GetPropNames():
        mol.SetProp(prop, mol.GetProp(prop))
    writer.write(mol)
writer.close()

# Convert each molecule to PDBQT using its name
supplier = Chem.SDMolSupplier("mols_with_H.sdf", removeHs=False)

for i, mol in enumerate(supplier):
    if mol is None:
        continue
    name = mol.GetProp("_Name").strip().replace(" ", "_")
    sdf_filename = f"docking_files/{name}.sdf"
    pdbqt_filename = f"docking_files/{name}.pdbqt"

    # Write individual SDF
    writer = Chem.SDWriter(sdf_filename)
    writer.write(mol)
    writer.close()

    # Convert to PDBQT with Open Babel
    # xh tells openbabel to keep hydrogens
    os.system(f"obabel {sdf_filename} -O {pdbqt_filename} -xh")

print("Docking preparation complete. Files are ready in 'docking_files' directory.")
