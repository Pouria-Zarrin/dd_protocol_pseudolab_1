#!/usr/bin/env python3
"""
This script runs all bash executable scripts in a specified directory.
Usage:
    python run_all_bash_scripts.py --path_to_scripts <directory_with_bash_scripts>
"""
import argparse
import os
import subprocess
import sys
import concurrent.futures


def is_bash_executable(filepath):
    """Check if a file is a bash script and executable."""
    if not os.path.isfile(filepath):
        return False
    if not os.access(filepath, os.X_OK):
        return False
    with open(filepath, 'r') as f:
        first_line = f.readline()
        return first_line.startswith('#!') and 'bash' in first_line


def run_script(script):
    try:
        result = subprocess.run([script], capture_output=True, text=True, check=True)
        return (script, True, result.stdout, result.stderr)
    except subprocess.CalledProcessError as e:
        return (script, False, e.stdout, e.stderr)
    except Exception as ex:
        return (script, False, '', str(ex))


def replace_progressive_docking_path(directory, progressive_docking_path):
    """Replace 'python progressive_docking.py' with 'python <progressive_docking_path>' in all .sh files in the directory."""
    for fname in os.listdir(directory):
        if fname.endswith('.sh'):
            fpath = os.path.join(directory, fname)
            with open(fpath, 'r') as file:
                content = file.read()
            # Replace only the command, not other possible mentions
            new_content = content.replace('python -u progressive_docking.py', f'python -u {progressive_docking_path}')
            new_content = new_content.replace('python progressive_docking.py', f'python {progressive_docking_path}')
            if new_content != content:
                with open(fpath, 'w') as file:
                    file.write(new_content)


def run_bash_scripts_in_dir(directory, progressive_docking_path, max_workers=4):
    # First, update all .sh files to use the correct path
    replace_progressive_docking_path(directory, progressive_docking_path)
    scripts = []
    for fname in sorted(os.listdir(directory)):
        fpath = os.path.join(directory, fname)
        if is_bash_executable(fpath):
            scripts.append(fpath)
    if not scripts:
        print(f"No bash executable scripts found in {directory}")
        return
    print(f"Found {len(scripts)} bash scripts to execute in {directory}:")
    for i, script in enumerate(scripts, 1):
        print(f"  {i}. {os.path.basename(script)}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_script, script): script for script in scripts}
        for i, future in enumerate(concurrent.futures.as_completed(futures), 1):
            script, success, stdout, stderr = future.result()
            print(f"\n[{i}/{len(scripts)}] Finished: {script}")
            if success:
                print(f"SUCCESS: {script}")
                print("--- STDOUT ---\n" + (stdout or ''))
                if stderr:
                    print("--- STDERR ---\n" + stderr)
            else:
                print(f"ERROR: {script}")
                print("--- STDOUT ---\n" + (stdout or ''))
                print("--- STDERR ---\n" + (stderr or ''))


def main():
    parser = argparse.ArgumentParser(
        description="Run all bash executable scripts in a directory with progress reporting."
    )
    parser.add_argument(
        "--path_to_scripts",
        required=True,
        help="Path to the directory containing bash scripts to execute."
    )
    parser.add_argument(
        "--max_workers",
        type=int,
        default=4,
        help="Maximum number of scripts to run in parallel (default: 4)."
    )
    parser.add_argument(
        "--progressive_docking_path",
        required=True,
        help="The path to replace in the bash scripts for progressive docking."
    )
    args = parser.parse_args()
    if not os.path.isdir(args.path_to_scripts):
        print(f"ERROR: {args.path_to_scripts} is not a valid directory.")
        sys.exit(1)
    run_bash_scripts_in_dir(args.path_to_scripts, args.progressive_docking_path, max_workers=args.max_workers)


if __name__ == "__main__":
    main()
