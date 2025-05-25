import os
import argparse
from openbabel import openbabel, pybel

def find_pdbqt_files(root_dir):
    pdbqt_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith('.pdbqt'):
                pdbqt_files.append(os.path.join(dirpath, filename))
    return pdbqt_files

def check_pdbqt_files(ligand_dir):
    pdbqt_files = find_pdbqt_files(ligand_dir)
    print(f"Found {len(pdbqt_files)} .pdbqt files to check.")
    corrupted_or_empty_files = []
    for pdbqt_file in pdbqt_files:
        # Check if file is empty
        if os.path.getsize(pdbqt_file) == 0:
            print(f"[EMPTY FILE] {pdbqt_file}")
            corrupted_or_empty_files.append(pdbqt_file)
            continue
        try:
            # Try to read the file with Open Babel
            mols = list(pybel.readfile('pdbqt', pdbqt_file))
            if not mols:
                print(f"[OPENBABEL EMPTY] {pdbqt_file}")
                corrupted_or_empty_files.append(pdbqt_file)
        except Exception as e:
            print(f"[OPENBABEL ERROR] {pdbqt_file}: {e}")
            corrupted_or_empty_files.append(pdbqt_file)
    print("\nCorrupted, unreadable, or empty files:")
    for f in corrupted_or_empty_files:
        print(f)
    print(f"Total problematic files: {len(corrupted_or_empty_files)}")

def main():
    parser = argparse.ArgumentParser(description="Check .pdbqt files for parse errors using Open Babel.")
    parser.add_argument('--ligands_pdbqt_dir', required=True, help='Path to the directory containing .pdbqt ligand files (may have subdirectories).')
    args = parser.parse_args()
    check_pdbqt_files(args.ligands_pdbqt_dir)

if __name__ == "__main__":
    main()
