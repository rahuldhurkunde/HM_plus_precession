import h5py
import numpy as np

def write_zero_spins(filename):
	hf = h5py.File(filename, 'r+')
	nsignals = len(hf['mass1'])
	for n in range(len(hf['mass1'])):
		hf['spin1x'][n] = 0.0
		hf['spin1y'][n] = 0.0
		hf['spin2x'][n] = 0.0
		hf['spin2y'][n] = 0.0
	#hf.create_dataset('spin1x', data=np.zeros(np.shape(hf['mass1'])))	
	hf.close()

sets = 51
for k in range(sets):
	filename = '%s.hdf' %k
	
	write_zero_spins(filename)
	#f = h5py.File(filename, 'r')
	#print(np.shape(f['spin1z']), f['spin1z'][0])	
