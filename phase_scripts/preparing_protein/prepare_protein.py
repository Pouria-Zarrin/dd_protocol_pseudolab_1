import os
import requests
import json
from MDAnalysis import Universe
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from openmm.app import Modeller
from openmm.app import Element
from openmm.app import ForceField
from openmm.app import PDBxFile
from openmm.app import PDBFile
from openmm.app import Modeller
from openmm.app import Element
from openmm.app import ForceField
from openmm.app import PDBxFile
from openmm.app import PDBFile
from openmm.app import Modeller
from openmm.app import Element
from openmm.app import ForceField
from openmm.app import PDBxFile
import subprocess


def download_pdb_by_id(pdb_id, out_path):
    """
    Download the PDB structure for a given PDB ID from RCSB PDB.
    """
    pdb_request = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
    pdb_request.raise_for_status()
    with open(out_path, 'w') as f:
        f.write(pdb_request.text)
    return pdb_id


def split_pdb_residues(pdb_path, output_dir, ligand_name):
    """
    Split the PDB into protein, ligand, and water using MDAnalysis and save as separate files.
    """
    u = Universe(pdb_path)
    protein = u.select_atoms("protein")
    ligand = u.select_atoms(f"resname {ligand_name}")
    water = u.select_atoms("resname HOH")
    protein.write(os.path.join(output_dir, "protein.pdb"))
    ligand.write(os.path.join(output_dir, "ligand.pdb"))
    water.write(os.path.join(output_dir, "water.pdb"))


def fix_protein_with_pdbfixer(input_pdb, output_dir):
    """
    Use PDBFixer to fix the protein structure and save the fixed PDB and PQR files.
    """
    os.makedirs(output_dir, exist_ok=True)
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()
    fixed_pdb = os.path.join(output_dir, "fixed_protein.pdb")
    with open(fixed_pdb, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    # Save as PQR (requires openmmforcefields or pdb2pqr, not directly in PDBFixer)
    # Here, we just save the fixed PDB. For PQR, you may need to use external tools.
    return fixed_pdb


def convert_pdb_to_pdbqt(pdb_path, pdbqt_dir):
    """
    Convert a PDB file to PDBQT using Open Babel.
    """
    os.makedirs(pdbqt_dir, exist_ok=True)
    pdbqt_path = os.path.join(pdbqt_dir, os.path.splitext(os.path.basename(pdb_path))[0] + ".pdbqt")
    subprocess.run(["obabel", pdb_path, "-O", pdbqt_path], check=True)
    return pdbqt_path


def main(pdb_id, ligand_name):
    base_dir = os.path.dirname(os.path.abspath(__file__))
    pdb_path = os.path.join(base_dir, "original_protein_from_pdb.pdb")
    print(f"Downloading PDB for ID {pdb_id}...")
    pdb_id = download_pdb_by_id(pdb_id, pdb_path)
    print(f"Downloaded {pdb_id} to {pdb_path}")

    split_dir = os.path.join(base_dir, "split_residues")
    os.makedirs(split_dir, exist_ok=True)
    print("Splitting PDB into protein, ligand, and water...")
    split_pdb_residues(pdb_path, split_dir, ligand_name)

    fixed_dir = os.path.join(base_dir, "fixed protein")
    print("Fixing protein with PDBFixer...")
    fixed_pdb = fix_protein_with_pdbfixer(os.path.join(split_dir, "protein.pdb"), fixed_dir)

    pdbqt_dir = os.path.join(base_dir, "pdbqt")
    print("Converting fixed PDB to PDBQT with Open Babel...")
    convert_pdb_to_pdbqt(fixed_pdb, pdbqt_dir)
    print("All steps complete.")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python prepare_protein.py <PDB_ID> <LIGAND_RESNAME>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
