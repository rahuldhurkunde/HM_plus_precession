import numpy as np
import h5py
import sys
import argparse as ap

parser = ap.ArgumentParser()
parser.add_argument('--tb_path',
        required = True,
        help ='specify TB path (default file name -- combined_bank.hdf)')
args, remaining_args = parser.parse_known_args()

tb_file = '%s/combined_bank.hdf' %args.tb_path
print(" \t \t Bank given", tb_file)
f = h5py.File(tb_file, 'r+')

tau0_temp = np.array(f['tau0'][:])
indices = np.argsort(tau0_temp)

mass1 = f['mass1'][:]
mass2 = f['mass2'][:]
spin1z = f['spin1z'][:]
spin2z = f['spin2z'][:]
tau0 = f['tau0'][:]
tau3 = f['tau3'][:]

newfile = '%s/sorted_bank.hdf' %args.tb_path
with h5py.File(newfile, 'w') as hf:
	hf.create_dataset("mass1", data=mass1[indices])
	hf.create_dataset("mass2", data=mass2[indices])
	hf.create_dataset("spin1z", data=spin1z[indices])
	hf.create_dataset("spin2z", data=spin2z[indices])
	hf.create_dataset("tau0", data=tau0[indices])
	hf.create_dataset("tau3", data=tau3[indices])
hf.close()
print('\t Sorting complete and stored in', newfile)
