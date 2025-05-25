import os
import argparse
from vina import Vina
import subprocess

# Set up paths (edit these as needed)
LIGAND_PATH = "projects/project_0/iteration_1/smile/docking_ready/docking_ready_pdbqt_train/ZINC000607560302.pdbqt"
RECEPTOR_PATH = "phase_scripts/preparing_protein/pdbqt/fixed_protein.pdbqt"
OUTPUT_DIR = "projects/project_0/iteration_1/smile/docking_results/single_ligand_test"
OUT_PDBQT = os.path.join(OUTPUT_DIR, "ZINC000607560302_docked.pdbqt")
OUT_SDF = os.path.join(OUTPUT_DIR, "ZINC000607560302_docked.sdf")

# Box and center (edit as needed)
CENTER = [-24.62200164794922, -0.3884687125682831, -10.92924976348877]
BOX = [21.163999557495117, 16.291000366210938, 19.618000030517578]
NUM_POSES = 5

os.makedirs(OUTPUT_DIR, exist_ok=True)

v = Vina(sf_name='vina')
v.set_receptor(RECEPTOR_PATH)
v.compute_vina_maps(center=CENTER, box_size=BOX)
v.set_ligand_from_file(LIGAND_PATH)
v.dock()
v.write_poses(OUT_PDBQT, n_poses=NUM_POSES, overwrite=True)

# Convert to SDF using Open Babel
result = subprocess.run([
    'obabel', OUT_PDBQT, '-O', OUT_SDF
], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if result.returncode == 0:
    print(f"SDF result written: {OUT_SDF}")
else:
    print(f"Open Babel failed to convert {OUT_PDBQT} to SDF. Error: {result.stderr.decode()}")
