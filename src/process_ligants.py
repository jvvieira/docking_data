import os
from xml.parsers.expat import model
import pandas as pd
import numpy as np
from scipy import io
from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select

parser = PDBParser(QUIET=True)

files = os.listdir("./filtered_data")
print(f"Files in /filtered_data: {len(files)}")

def get_ligands_from_pdb(pdb_file):
    # 1. Parse the structure
    structure = parser.get_structure(pdb_file, f"{pdb_file}")
  
    # 2. Separate protein and ligand/ligands
    protein_atoms = []
    ligand_atoms = []
    for atom in structure.get_atoms():
        if(atom.parent.parent.get_id() == 'A'): # Example: only consider chain A as protein
            protein_atoms.append(atom)
        elif(atom.parent.parent.get_id() == 'B'): # Example: only consider chain B as ligand
            ligand_atoms.append(atom)
            
    # print(f"Protein atoms: {len(protein_atoms)}, Ligand atoms: {len(ligand_atoms)}")

    # 3. Find intersection (atoms within 5 Angstroms)
    ns = NeighborSearch(protein_atoms)
    nearby_residues = set()
    for ligand_atom in ligand_atoms:
        # Find protein residues within 5A of this ligand atom
        for neighbor in ns.search(ligand_atom.coord, 5.0, level='A'):
            nearby_residues.add(neighbor)

    return nearby_residues


ligant1 = get_ligands_from_pdb("./filtered_data/RXR_TS15_model4_prep_0028.pdb")
print(len(ligant1))
# Iterate through found atoms
for atom in ligant1:
    print(f"Residue: {atom.name} {atom}")