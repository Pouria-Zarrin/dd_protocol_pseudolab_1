import os
import glob
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops
import subprocess

INPUT_DIR = "projects/project_0/iteration_1/smile"
OUTPUT_SUFFIX = "_prepared_ligand.smi"


def optimize_smiles_file(input_path, output_path, pdbqt_dir):
    with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            smiles, name = parts[0], parts[1]
            # Extract ZINC ID (assume it's the part before the first underscore)
            zinc_id = name.split('_')[0] if '_' in name else name
            mol = Chem.MolFromSmiles(smiles, sanitize=True)
            if mol is None:
                print(f"Skipping molecule with invalid SMILES: {name} (ZINC ID: {zinc_id})")
                continue
            mol = Chem.AddHs(mol)
            # Try constrained embedding and geometry optimization
            try:
                AllChem.EmbedMolecule(mol, useRandomCoords=True)
                AllChem.UFFOptimizeMolecule(mol)
                smiles_opt = Chem.MolToSmiles(mol, isomericSmiles=True)
                outfile.write(f"{smiles_opt} {name}\n")
                # Save as MOL file for Open Babel conversion
                mol_file = output_path.replace('.smi', f'_{zinc_id}.mol')
                Chem.MolToMolFile(mol, mol_file)
                # Convert MOL to PDBQT using Open Babel with timeout
                pdbqt_file = os.path.join(pdbqt_dir, f"{zinc_id}.pdbqt")
                try:
                    result = subprocess.run(
                        ["obabel", mol_file, "-O", pdbqt_file, "--gen3d"],
                        stderr=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        timeout=30  # seconds
                    )
                    if result.returncode != 0 or not os.path.exists(pdbqt_file):
                        with open("obabel_error.log", "a") as errlog:
                            errlog.write(f"Open Babel failed for {name} (ZINC ID: {zinc_id}):\n{result.stderr.decode()}\n")
                        print(f"FAILED: Open Babel failed for {name} (ZINC ID: {zinc_id}). See obabel_error.log for details.")
                    else:
                        print(f"SUCCESS: Open Babel conversion succeeded for {name} (ZINC ID: {zinc_id})")
                except subprocess.TimeoutExpired:
                    with open("obabel_error.log", "a") as errlog:
                        errlog.write(f"Open Babel timed out for {name} (ZINC ID: {zinc_id})\n")
                    print(f"TIMEOUT: Open Babel timed out for {name} (ZINC ID: {zinc_id}). Skipping.")
                finally:
                    if os.path.exists(mol_file):
                        os.remove(mol_file)
            except Exception as e:
                # If optimization fails, write the original SMILES
                outfile.write(f"{smiles} {name}\n")
                print(f"ERROR: RDKit optimization failed for {name} (ZINC ID: {zinc_id}): {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Optimize SMILES using RDKit and save to output directory."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input .smi file or directory containing .smi files",
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to save prepared ligand .smi files",
    )
    args = parser.parse_args()

    input_path = args.input
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    if os.path.isdir(input_path):
        smi_files = glob.glob(os.path.join(input_path, "*.smi"))
    else:
        smi_files = [input_path]

    for smi_file in smi_files:
        base = os.path.splitext(os.path.basename(smi_file))[0]
        output_file = os.path.join(output_dir, base + OUTPUT_SUFFIX)
        # Determine set type (train, test, valid) from file name
        set_type = base.split('_')[0]
        pdbqt_dir = os.path.join(output_dir, f"docking_ready_pdbqt_{set_type}")
        os.makedirs(pdbqt_dir, exist_ok=True)
        optimize_smiles_file(smi_file, output_file, pdbqt_dir)
        print(f"Processed: {smi_file} -> {output_file} and PDBQT files in {pdbqt_dir}")


if __name__ == "__main__":
    main()
