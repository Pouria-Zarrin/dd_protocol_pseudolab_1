import glob
import time
import numpy as np
import pandas as pd
import pickle
from contextlib import closing
from multiprocessing import Pool
import multiprocessing
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import Chem
from functools import partial
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-sfp','--smile_folder_path',help='name of the folder with prepared smiles',required=True)
parser.add_argument('-fn','--folder_name',help='name of morgan fingerprint folder',required=True)
parser.add_argument('-tp','--tot_process',help='number of cores',required=True)

io_args = parser.parse_args()
sfp = io_args.smile_folder_path
fn = io_args.folder_name
t_pos = int(io_args.tot_process)

# Check directory accessibility
if not os.path.exists(sfp):
    print(f"ERROR: The smiles folder path '{sfp}' does not exist or is not accessible.")
    exit(1)
else:
    print(f"Smiles folder path '{sfp}' is accessible. Files: {os.listdir(sfp)}")

if not os.path.exists(fn):
    print(f"NOTE: The output folder '{fn}' does not exist yet. It will be created.")
else:
    print(f"Output folder '{fn}' is accessible. Files: {os.listdir(fn)}")

def morgan_fingp(fname):
    nbits=1024
    radius=2
    fsplit = fname.split('/')[-1]
    ref2  = open(fn+'/'+fsplit,'a')
    with open(fname,'r') as ref:
        for line in ref:
            try:
                smile,zin_id = line.rstrip().split()
                arg = np.zeros((nbits,), dtype=int)
                DataStructs.ConvertToNumpyArray(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smile),radius,nBits=nbits,useChirality=True),arg)
                ref2.write((',').join([zin_id]+[str(elem) for elem in np.where(arg==1)[0]]))
                ref2.write('\n')
            except Exception as e:
                print(f"Error processing line: {line.strip()} | Error: {e}")
                pass

files = []
for f in glob.glob(sfp+'/*.smi'):
    files.append(f)

try:
    os.mkdir(fn)
except:
    pass

t_f = len(files)
t = time.time()
with closing(Pool(np.min([multiprocessing.cpu_count(),t_pos]))) as pool:
    pool.map(morgan_fingp,files)
print(time.time()-t)
