import argparse
import os
import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor

# Argument parsing
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project', required=True, help='Project name (e.g., project_0)')
    parser.add_argument('--iteration', required=True, type=int, help='Iteration number (e.g., 1)')
    parser.add_argument('--protein', required=True, help='Protein name')
    parser.add_argument('--morgan_library', required=True, help='Morgan library directory')
    parser.add_argument('--file_path', required=True, help='File path for output')
    parser.add_argument('--python_file', required=True, help='Python file to run in the .sh scripts')
    parser.add_argument('--max_workers', type=int, default=4, help='Number of parallel jobs')
    return parser.parse_args()

def make_executable(path):
    os.chmod(path, 0o755)

def edit_sh_file(sh_path, python_file, fn, protein, morgan_library, file_path):
    with open(sh_path, 'r') as f:
        lines = f.readlines()
    new_lines = []
    for line in lines:
        if line.strip().startswith('python -u'):
            # Replace the command line with new arguments
            new_cmd = f"python -u {python_file} -fn {fn} -protein {protein} -it 1 -mdd {morgan_library} -file_path {file_path}\n"
            new_lines.append(new_cmd)
        else:
            new_lines.append(line)
    with open(sh_path, 'w') as f:
        f.writelines(new_lines)

def run_sh_file(sh_path):
    subprocess.run([sh_path], check=True)

def main():
    args = parse_args()
    sh_dir = f"projects/{args.project}/iteration_{args.iteration}/simple_job_predictions"
    sh_files = sorted(glob.glob(os.path.join(sh_dir, '*.sh')))
    if not sh_files:
        print(f"No .sh files found in {sh_dir}")
        return
    # Make all .sh files executable
    for sh in sh_files:
        make_executable(sh)
    # Edit each .sh file with the provided arguments
    for sh in sh_files:
        # Extract the .smi filename from the command in the .sh file
        with open(sh, 'r') as f:
            for line in f:
                if line.strip().startswith('python -u'):
                    parts = line.strip().split()
                    try:
                        fn_idx = parts.index('-fn') + 1
                        fn = parts[fn_idx]
                    except Exception:
                        fn = ''
                    break
        edit_sh_file(sh, args.python_file, fn, args.protein, args.morgan_library, args.file_path)
    # Run all .sh files in parallel
    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = [executor.submit(run_sh_file, sh) for sh in sh_files]
        for future in futures:
            try:
                future.result()
            except Exception as e:
                print(f"Error running script: {e}")

if __name__ == '__main__':
    main()
