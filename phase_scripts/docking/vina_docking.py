import os
import argparse
import numpy as np
import MDAnalysis as mda
from vina import Vina
import subprocess


def get_ligand_and_pocket_center(protein_path, ligand_resname, box_expansion):
    u = mda.Universe(protein_path)
    ligand = u.select_atoms(f'resname {ligand_resname}')
    protein = u.select_atoms('protein')
    if len(ligand) == 0:
        raise ValueError(f"No ligand with resname {ligand_resname} found in structure.")
    ligand_coords = ligand.positions
    min_corner = ligand_coords.min(axis=0) - box_expansion
    max_corner = ligand_coords.max(axis=0) + box_expansion
    center = ligand_coords.mean(axis=0)
    box = max_corner - min_corner
    return center.tolist(), box.tolist()


def dock_and_save(vina_obj, pdbqt_dir, output_dir, num_poses=5):
    for subdir in os.listdir(pdbqt_dir):
        subdir_path = os.path.join(pdbqt_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue
        docking_result_dir = os.path.join(output_dir, f"{subdir}_docking_result")
        os.makedirs(docking_result_dir, exist_ok=True)
        for file in os.listdir(subdir_path):
            if file.endswith('.pdbqt'):
                ligand_path = os.path.join(subdir_path, file)
                print(f"Docking ligand: {ligand_path}")
                vina_obj.set_ligand_from_file(ligand_path)
                vina_obj.dock()
                out_pdbqt = os.path.join(docking_result_dir, file)
                vina_obj.write_poses(out_pdbqt, n_poses=num_poses, overwrite=True)
                print(f"Docking result written: {out_pdbqt}")
                # Convert the best pose to SDF using Open Babel
                sdf_out = os.path.join(docking_result_dir, file.replace('.pdbqt', '.sdf'))
                try:
                    result = subprocess.run([
                        'obabel', out_pdbqt, '-O', sdf_out
                    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    if result.returncode == 0:
                        print(f"SDF result written: {sdf_out}")
                    else:
                        print(f"Open Babel failed to convert {out_pdbqt} to SDF. Error: {result.stderr.decode()}")
                except Exception as e:
                    print(f"Failed to convert {out_pdbqt} to SDF: {e}")


def main():
    parser = argparse.ArgumentParser(description="Automated Vina docking workflow.")
    parser.add_argument('--original_protein_path', required=True, help='Path to the original protein structure (PDB)')
    parser.add_argument('--ligand_resname', required=True, help='Residue name of the ligand in the structure')
    parser.add_argument('--box_expansion', type=float, default=4.0, help='Box expansion in angstroms')
    parser.add_argument('--fixed_protein_pdbqt', required=True, help='Path to the fixed protein PDBQT file')
    parser.add_argument('--ligands_pdbqt_dir', required=True, help='Directory containing subdirectories with ligand PDBQT files')
    parser.add_argument('--output_dir', required=True, help='Directory to save docking results')
    parser.add_argument('--exhaustiveness', type=int, default=8, help='Vina exhaustiveness')
    parser.add_argument('--num_poses', type=int, default=5, help='Number of docking poses to save')
    args = parser.parse_args()

    center, box = get_ligand_and_pocket_center(args.original_protein_path, args.ligand_resname, args.box_expansion)
    print(f"Pocket center: {center}")
    print(f"Box size: {box}")

    # Replace 'TITLE' and 'CRYST1' with 'REMARK' in the fixed protein PDBQT file before using it
    with open(args.fixed_protein_pdbqt, 'r') as f:
        file_content = f.readlines()
    # Remove lines starting with ROOT, ENDROOT, BRANCH, ENDBRANCH, TORSDOF
    cleaned_lines = [line for line in file_content if not (line.strip().startswith('ROOT') or line.strip().startswith('ENDROOT') or line.strip().startswith('BRANCH') or line.strip().startswith('ENDBRANCH') or line.strip().startswith('TORSDOF'))]
    # Replace TITLE and CRYST1 with REMARK
    cleaned_lines = [line.replace('TITLE', 'REMARK').replace('CRYST1', 'REMARK') for line in cleaned_lines]
    with open(args.fixed_protein_pdbqt, 'w') as f:
        f.writelines(cleaned_lines)

    v = Vina(sf_name='vina')
    v.set_receptor(args.fixed_protein_pdbqt)
    v.compute_vina_maps(center=center, box_size=box)
    dock_and_save(v, args.ligands_pdbqt_dir, args.output_dir, num_poses=args.num_poses)

if __name__ == "__main__":
    main()
