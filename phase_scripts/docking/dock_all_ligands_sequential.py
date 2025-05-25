import os
from vina import Vina
import subprocess

LIGAND_DIR = "projects/project_0/iteration_1/smile/docking_ready/docking_ready_pdbqt_train"
RECEPTOR_PATH = "phase_scripts/preparing_protein/pdbqt/fixed_protein.pdbqt"
OUTPUT_DIR = "projects/project_0/iteration_1/smile/docking_results/docking_ready_pdbqt_train_docking_result"
CENTER = [-24.62200164794922, -0.3884687125682831, -10.92924976348877]
BOX = [21.163999557495117, 16.291000366210938, 19.618000030517578]
NUM_POSES = 5

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Precompute already docked ligands
already_docked = set()
for file in os.listdir(OUTPUT_DIR):
    if file.endswith('.sdf') and os.path.getsize(os.path.join(OUTPUT_DIR, file)) > 0:
        already_docked.add(file.replace('.sdf', '.pdbqt'))

ligand_files = [f for f in os.listdir(LIGAND_DIR) if f.endswith('.pdbqt')]
for ligand_file in ligand_files:
    if ligand_file in already_docked:
        print(f"Skipping already docked ligand: {ligand_file}")
        continue
    ligand_path = os.path.join(LIGAND_DIR, ligand_file)
    out_pdbqt = os.path.join(OUTPUT_DIR, ligand_file.replace('.pdbqt', '_docked.pdbqt'))
    out_sdf = os.path.join(OUTPUT_DIR, ligand_file.replace('.pdbqt', '_docked.sdf'))
    try:
        v = Vina(sf_name='vina')
        v.set_receptor(RECEPTOR_PATH)
        v.compute_vina_maps(center=CENTER, box_size=BOX)
        v.set_ligand_from_file(ligand_path)
        v.dock()
        v.write_poses(out_pdbqt, n_poses=NUM_POSES, overwrite=True)
        result = subprocess.run([
            'obabel', out_pdbqt, '-O', out_sdf
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            print(f"Docked and converted: {ligand_file}")
        else:
            print(f"Docked but Open Babel failed for: {ligand_file}. Error: {result.stderr.decode()}")
    except Exception as e:
        print(f"Error docking {ligand_file}: {e}")
