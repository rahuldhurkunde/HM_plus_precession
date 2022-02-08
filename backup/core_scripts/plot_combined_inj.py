import h5py as hf
import numpy as np
import matplotlib.pyplot as plt
import functions as func
import pycbc
from pycbc import conversions


tb = func.read_tb("/work/rahul.dhurkunde/HM_and_precession/banks/parallel/aligned_bank/combined_bank.hdf", 30.0)
tb_m1 = [x.m1 for x in tb]
tb_m2 = [x.m2 for x in tb]

nsplits = 250 
for k in range(nsplits):
	filename = "50000_inj/nonaligned_injections/%s.hdf" %k
	f = hf.File(filename, 'r')
	m1 = np.array(f['mass1'])
	m2 = np.array(f['mass2'])
	spin1x = np.array(f['spin1x'])
	spin1y = np.array(f['spin1y'])
	spin2x = np.array(f['spin2x'])
	spin2y = np.array(f['spin2y'])
	inclination = np.degrees(np.array(f['inclination']))
	chi_p = conversions.chi_p(m1[:], m2[:], spin1x[:], spin1y[:], spin2x[:], spin2y[:])
	plt.scatter(chi_p, inclination)
print(len(m1))
#plt.plot(tb_m1, tb_m2, '.', label='templates', alpha = 0.1)
plt.show()
