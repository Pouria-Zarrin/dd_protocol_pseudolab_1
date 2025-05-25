import os
import time

morgan_dir = '../../morgan_library'

for fname in os.listdir(morgan_dir):
    fpath = os.path.join(morgan_dir, fname)
    if os.path.isfile(fpath):
        mod_time = os.path.getmtime(fpath)
        print(f"{fname}: {time.ctime(mod_time)}")
