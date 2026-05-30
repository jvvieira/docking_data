import os
from xml.parsers.expat import model
import pandas as pd
import numpy as np
from scipy import io
from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Select

parser = PDBParser(QUIET=True)

files = os.listdir("./filtered_data")

def get_ligands_from_pdb(pdb_file):
    # 1. Parse the structure
    structure = parser.get_structure(pdb_file, f"{pdb_file}")
  
    # 2. Separate protein and ligand/ligands
    protein_residues = []
    ligand_atoms = []
    for residue in structure.get_residues():
        if(residue.parent.get_id() == 'A'): # Example: only consider chain A as protein
            protein_residues.append(residue)
        elif(residue.parent.get_id() == 'B'): # Example: only consider chain B as ligand
            ligand_atoms.append(residue)
            
    # print(f"Protein residues: {len(protein_residues)}, Ligand residues: {len(ligand_atoms)}")

    # 3. Find intersection (residues within 5 Angstroms)
    ns = NeighborSearch([atom for res in protein_residues for atom in res.get_atoms()])
    nearby_residues = set()
    for ligand in ligand_atoms:
        for atom in ligand.get_atoms():
            neighbors = ns.search(atom.coord, 5.0)  # Search for neighbors within 5 Angstroms
            for neighbor in neighbors:
                nearby_residues.add(neighbor.get_parent())  # Add the residue of the neighboring atom
                
    return nearby_residues

files = os.listdir("./filtered_data")
size = len(files)

finaldata = pd.DataFrame(columns=['protein', 'ligand', 'model', 'pose', 'residue', 'description'])

for file in files:
    # print(f"Processing file: {file} ({files.index(file)+1}/{size})")
    ligand = get_ligands_from_pdb(f"./filtered_data/{file}")
    for residue in ligand:
        new_data = {
            'protein': file.split("_")[0],
            'ligand': file.split("_")[1],
            'model': file.split("_")[2],
            'pose': file.split("_")[4].replace(".pdb", ""),
            'description': file.replace(".pdb", ""),
            'residue': residue.get_resname()
        }
        finaldata = pd.concat([finaldata, pd.DataFrame([new_data])], ignore_index=True)

finaldata.to_csv('./data/ligand_residues.csv', index=False)