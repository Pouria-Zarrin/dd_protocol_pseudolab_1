import os
import argparse
import numpy as np
import MDAnalysis as mda
from vina import Vina
import subprocess
import concurrent.futures

# Get the directory where this script is located
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ERROR_LOG_PATH = os.path.join(SCRIPT_DIR, 'error_log.txt')


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


def dock_single_ligand(vina_obj_args, ligand_path, out_pdbqt, num_poses, sdf_out):
    try:
        from vina import Vina
        import subprocess
        print(f"[DOCKING] Processing ligand: {ligand_path}")
        with open(ERROR_LOG_PATH, 'a') as logf:
            logf.write(f"[INFO] Starting docking for ligand: {ligand_path}\n")
        v = Vina(sf_name='vina')
        try:
            v.set_receptor(vina_obj_args['receptor'])
        except BaseException as e:
            error_msg = f"Failed to set receptor for {ligand_path}: {e}"
            print(f"[ERROR] {error_msg}")
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[Docking Error] Ligand: {ligand_path}, Step: set_receptor, Exception: {e}\n")
            return
        try:
            v.compute_vina_maps(center=vina_obj_args['center'], box_size=vina_obj_args['box'])
        except BaseException as e:
            error_msg = f"Failed to compute vina maps for {ligand_path}: {e}"
            print(f"[ERROR] {error_msg}")
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[Docking Error] Ligand: {ligand_path}, Step: compute_vina_maps, Exception: {e}\n")
            return
        print(f"Docking ligand: {ligand_path}")
        try:
            v.set_ligand_from_file(ligand_path)
        except BaseException as e:
            error_msg = f"Failed to set ligand from file for {ligand_path}: {e}"
            print(f"[ERROR] {error_msg}")
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[Docking Error] Ligand: {ligand_path}, Step: set_ligand_from_file, Exception: {e}\n")
            return
        try:
            v.dock()
        except BaseException as e:
            error_msg = f"Docking failed for {ligand_path}: {e}"
            print(f"[ERROR] {error_msg}")
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[Docking Error] Ligand: {ligand_path}, Step: dock, Exception: {e}\n")
            return
        try:
            v.write_poses(out_pdbqt, n_poses=num_poses, overwrite=True)
            print(f"Docking result written: {out_pdbqt}")
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[INFO] Docking completed for ligand: {ligand_path}, output: {out_pdbqt}\n")
        except BaseException as e:
            error_msg = f"Failed to write poses for {ligand_path}: {e}"
            print(f"[ERROR] {error_msg}")
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[Docking Error] Ligand: {ligand_path}, Step: write_poses, Exception: {e}\n")
            return
        try:
            result = subprocess.run([
                'obabel', out_pdbqt, '-O', sdf_out
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode == 0:
                print(f"SDF result written: {sdf_out}")
                with open(ERROR_LOG_PATH, 'a') as logf:
                    logf.write(f"[INFO] SDF conversion successful for ligand: {ligand_path}, sdf: {sdf_out}\n")
            else:
                error_msg = f"Open Babel failed to convert {out_pdbqt} to SDF. Error: {result.stderr.decode()}"
                print(error_msg)
                with open(ERROR_LOG_PATH, 'a') as logf:
                    logf.write(f"[OpenBabel Error] Ligand: {ligand_path}, Out: {out_pdbqt}, SDF: {sdf_out}, Error: {result.stderr.decode()}\n")
        except BaseException as e:
            error_msg = f"Failed to convert {out_pdbqt} to SDF: {e}"
            print(error_msg)
            with open(ERROR_LOG_PATH, 'a') as logf:
                logf.write(f"[OpenBabel Exception] Ligand: {ligand_path}, Out: {out_pdbqt}, SDF: {sdf_out}, Exception: {e}\n")
    except BaseException as e:
        error_msg = f"Docking failed for {ligand_path}: {e}"
        print(f"[ERROR] {error_msg}")
        with open(ERROR_LOG_PATH, 'a') as logf:
            logf.write(f"[Docking Error] Ligand: {ligand_path}, Out: {out_pdbqt}, SDF: {sdf_out}, Exception: {e}\n")


def dock_and_save(receptor, center, box, pdbqt_dir, output_dir, num_poses=5, max_workers=4):
    vina_obj_args = {
        'receptor': receptor,
        'center': center,
        'box': box
    }
    # Read failed ligands from previous runs
    failed_ligands = set()
    if os.path.exists('vina_docking_error.log'):
        with open('vina_docking_error.log', 'r') as logf:
            for line in logf:
                if ':' in line:
                    failed_ligands.add(line.split(':')[0].strip())
    tasks = []
    with open(ERROR_LOG_PATH, 'a') as logf:
        logf.write(f"[INFO] Starting batch docking. pdbqt_dir: {pdbqt_dir}, output_dir: {output_dir}\n")
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        for subdir in os.listdir(pdbqt_dir):
            subdir_path = os.path.join(pdbqt_dir, subdir)
            if not os.path.isdir(subdir_path):
                continue
            docking_result_dir = os.path.join(output_dir, f"{subdir}_docking_result")
            os.makedirs(docking_result_dir, exist_ok=True)
            for file in os.listdir(subdir_path):
                if file.endswith('.pdbqt'):
                    ligand_path = os.path.join(subdir_path, file)
                    out_pdbqt = os.path.join(docking_result_dir, file)
                    sdf_out = os.path.join(docking_result_dir, file.replace('.pdbqt', '.sdf'))
                    # Skip if SDF already exists and is non-empty
                    if os.path.exists(sdf_out) and os.path.getsize(sdf_out) > 0:
                        print(f"Skipping already docked ligand: {ligand_path}")
                        with open('error_log.txt', 'a') as logf:
                            logf.write(f"[INFO] Skipping already docked ligand: {ligand_path}\n")
                        continue
                    # Skip if ligand previously failed
                    if ligand_path in failed_ligands:
                        print(f"Skipping previously failed ligand: {ligand_path}")
                        with open('error_log.txt', 'a') as logf:
                            logf.write(f"[INFO] Skipping previously failed ligand: {ligand_path}\n")
                        continue
                    # Skip if ligand file is empty
                    if os.path.getsize(ligand_path) == 0:
                        print(f"Skipping empty ligand file: {ligand_path}")
                        with open('error_log.txt', 'a') as logf:
                            logf.write(f"[INFO] Skipping empty ligand file: {ligand_path}\n")
                        continue
                    print(f"Submitting ligand for docking: {ligand_path}")
                    with open('error_log.txt', 'a') as logf:
                        logf.write(f"[INFO] Submitting ligand for docking: {ligand_path}\n")
                    tasks.append(executor.submit(
                        dock_single_ligand, vina_obj_args, ligand_path, out_pdbqt, num_poses, sdf_out
                    ))
        for future in concurrent.futures.as_completed(tasks):
            try:
                future.result()
            except Exception as e:
                print(f"A docking task failed but continuing: {e}")
                with open('error_log.txt', 'a') as logf:
                    logf.write(f"[Batch Error] Exception in docking task: {e}\n")


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
    parser.add_argument('--max_workers', type=int, default=4, help='Maximum number of parallel docking processes (default: 4)')
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
    dock_and_save(
        args.fixed_protein_pdbqt,
        center,
        box,
        args.ligands_pdbqt_dir,
        args.output_dir,
        num_poses=args.num_poses,
        max_workers=args.max_workers
    )

if __name__ == "__main__":
    main()
