import os
import time

# Path to the output directory
morgan_dir = '../../morgan_library'

# Get initial set of files
if not os.path.exists(morgan_dir):
    print(f"Directory {morgan_dir} does not exist.")
    exit(1)

seen_files = set(os.listdir(morgan_dir))
print(f"Monitoring {morgan_dir} for new files...")

try:
    while True:
        current_files = set(os.listdir(morgan_dir))
        new_files = current_files - seen_files
        for f in new_files:
            print(f"New file written: {f}")
        seen_files = current_files
        time.sleep(2)
except KeyboardInterrupt:
    print("Stopped monitoring.")
