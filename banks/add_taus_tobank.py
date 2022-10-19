import h5py
import numpy as np
import argparse as ap
import pycbc
from pycbc import conversions

parser = ap.ArgumentParser()
parser.add_argument('--tb',
        required = True,
        help ='specify TB file')
parser.add_argument('--f_min',
        required = True, type=float,
        help ='specify lower frequency for tau computation')
args, remaining_args = parser.parse_known_args()
f_min = args.f_min

tb_file = '%s' %args.tb
hf = h5py.File(tb_file, 'r+')

mass1 = hf['mass1']
mass2 = hf['mass2']

print('TB --', args.tb, 'Total templates', len(mass1), 'f_min', f_min)
tau0 = []
tau3 = []
for n in range(len(mass1)):
	temp_tau0 = conversions.tau0_from_mass1_mass2(mass1[n], mass2[n], f_min)
	temp_tau3 = conversions.tau3_from_mass1_mass2(mass1[n], mass2[n], f_min)
	tau0.append(temp_tau0)
	tau3.append(temp_tau3)
	if (n % 1000 == 0):
		print(n)

with hf:
	hf.create_dataset("tau0", data=tau0)
	hf.create_dataset("tau3", data=tau3)
hf.close()

print('Added tau parameters into the TB')
