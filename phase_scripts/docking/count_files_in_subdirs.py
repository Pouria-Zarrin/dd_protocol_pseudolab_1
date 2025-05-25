import os

def count_files_in_subdirs(root_dir):
    print(f"Counting files in subdirectories of: {root_dir}")
    for subdir, dirs, files in os.walk(root_dir):
        if subdir == root_dir:
            continue  # Skip the root itself, only count subdirectories
        file_count = len([f for f in files if not f.startswith('.')])
        print(f"{subdir}: {file_count} files")

if __name__ == "__main__":
    dirs_to_check = [
        "projects/project_0/iteration_1/smile/docking_ready",
        "projects/project_0/iteration_1/smile/docking_results"
    ]
    for d in dirs_to_check:
        count_files_in_subdirs(d)
