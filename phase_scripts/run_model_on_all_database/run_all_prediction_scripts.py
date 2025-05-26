import os
import subprocess
import concurrent.futures

def make_executable(filepath):
    try:
        os.chmod(filepath, os.stat(filepath).st_mode | 0o111)
        print(f"Made executable: {filepath}")
    except Exception as e:
        print(f"Error making {filepath} executable: {e}")

def run_script(filepath):
    try:
        result = subprocess.run([filepath], capture_output=True, text=True, check=True)
        print(f"SUCCESS: {filepath}\n{result.stdout}")
        if result.stderr:
            print(f"STDERR from {filepath}:\n{result.stderr}")
        return (filepath, None)
    except subprocess.CalledProcessError as e:
        print(f"ERROR in {filepath} (exit code {e.returncode}):\nSTDOUT:\n{e.stdout}\nSTDERR:\n{e.stderr}")
        return (filepath, e.stderr)
    except Exception as e:
        print(f"UNEXPECTED ERROR in {filepath}: {e}")
        return (filepath, str(e))

def main():
    job_dir = os.path.join('projects', 'project_0', 'iteration_1', 'simple_job_predictions')
    if not os.path.isdir(job_dir):
        print(f"Directory not found: {job_dir}")
        return
    sh_files = [os.path.join(job_dir, f) for f in os.listdir(job_dir) if f.endswith('.sh')]
    if not sh_files:
        print("No .sh files found.")
        return
    for sh in sh_files:
        make_executable(sh)
    print(f"Running {len(sh_files)} scripts in parallel...")
    errors = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = {executor.submit(run_script, sh): sh for sh in sh_files}
        for future in concurrent.futures.as_completed(futures):
            filepath, error = future.result()
            if error:
                errors.append((filepath, error))
    if errors:
        print("\nErrors encountered:")
        for filepath, error in errors:
            print(f"{filepath}: {error}")
    else:
        print("All scripts executed successfully.")

if __name__ == "__main__":
    main()
