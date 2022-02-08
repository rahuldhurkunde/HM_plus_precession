import pycbc
from pycbc import waveform, psd, filter, conversions
import matplotlib.pyplot as plt
import numpy as np
import h5py
import injections as inj
import functions as func
import tqdm
import lal
import lalsimulation as lalsim
import time
lal.ClobberDebugLevel(0)

def load_injections(inj_dir, nsets, nsignals, f_min):
	sg = []
	for k in range(nsets):
		filename = '%s/%s.hdf' %(inj_dir, k)
		print('loaded', filename)
		sg = sg + inj.read_injections_HDF(filename, nsignals, f_min)

		m1 = [x.m1 for x in sg]
		m2 = [x.m2 for x in sg]
		spin1x = [x.s1x for x in sg]
		spin1y = [x.s1y for x in sg]
		spin2x = [x.s2x for x in sg]
		spin2y = [x.s2y for x in sg]
		inclination = np.degrees([x.inc for x in sg])
		chi_p = conversions.chi_p(m1[:], m2[:], spin1x[:], spin1y[:], spin2x[:], spin2y[:])
	return sg

def write_overlaps(filename, dot_products, sg_m1, sg_m2, mtotal, q):
	with h5py.File(filename, 'w') as f:
		f.create_dataset("overlap", data=dot_products)
		f.create_dataset("m1", data=sg_m1)
		f.create_dataset("m2", data=sg_m2)
		f.create_dataset("mtotal", data=mtotal)
		f.create_dataset("q", data=q)
	f.close()

def read_overlaps(filename):
	f = h5py.File(filename, 'r')
	overlap = np.array(f['overlap'])
	m1 = np.array(f['m1'])
	m2 = np.array(f['m2'])
	mtotal = np.array(f['mtotal'])
	q = np.array(f['q'])
	return overlap, m1, m2, mtotal, q

def get_thetaJN(mass1, mass2,
				spin1x, spin1y, spin1z,
				spin2x, spin2y, spin2z,
				inclination, fref, pref):

	thetaJN = []
	for k in range(len(mass1)): 
		temp = lalsim.SimIMRPhenomXPCalculateModelParametersFromSourceFrame(mass1[k], mass2[k], fref, pref, 
																					inclination[k], spin1x[k], spin1y[k], spin1z[k], 
																					spin2x[k], spin2y[k], spin2z[k], lalParams)
		thetaJN.append(temp[3])
	return np.array(thetaJN)
		

lalParams = lal.CreateDict()
f_min = 30.0
fref = 100.0
pref = 0.0
delta_f = 1.0/32
delta_t = 1.0/4096
approximant = 'IMRPhenomXPHM'

inj_dir = '/work/rahul.dhurkunde/HM_and_precession/injections/50000_inj/nonaligned_injections'
nsets = 250 
nsignals = 200

length = int(1.0/(delta_f*delta_t*2)) + 1
PSD = psd.analytical.aLIGOZeroDetHighPower(length, delta_f, f_min)

start = time.time()
sg = load_injections(inj_dir, nsets, nsignals, f_min)
end = time.time()
print("Inj loading time", end-start, 'Total injections', len(sg))

sg_m1 = [x.m1 for x in sg]
sg_m2 = [x.m2 for x in sg]
mtotal = sg_m1 + sg_m2
q = np.divide(sg_m1, sg_m2)
spin1x = [x.s1x for x in sg]
spin1y = [x.s1y for x in sg]
spin1z = [x.s1z for x in sg]
spin2x = [x.s2x for x in sg]
spin2y = [x.s2y for x in sg]
spin2z = [x.s2z for x in sg]
inclination = [x.inc for x in sg]
chi_p = conversions.chi_p(sg_m1[:], sg_m2[:], spin1x[:], spin1y[:], spin2x[:], spin2y[:])

thetaJN = get_thetaJN(sg_m1, sg_m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, inclination, fref, pref)

#dot_products = func.real_imag_dot_product(sg, PSD, f_min, delta_f, approximant, nsignals)
#write_overlaps("overlap.hdf", dot_products, sg_m1, sg_m2, mtotal, q)

overlap, m1, m2, mtotal, q = read_overlaps("overlap.hdf")
print(len(overlap), len(chi_p))

fig3, ax3 = plt.subplots()
#im3 = plt.scatter(m1, m2, c = overlap)
plt.xlabel('$\chi_{p}$')
plt.ylabel('inclination (deg)')
im3 = plt.scatter(chi_p, thetaJN, c = overlap)
fig3.colorbar(im3)
#plt.savefig('hplus_hcross.png', dpi = 600)
plt.show()


#indices = np.where(overlap > 0.)[0]
#m1_red = [m1[x] for x in indices]
#m2_red = [m2[x] for x in indices]
#mtotal_red = [mtotal[x] for x in indices]
#q_red = [q[x] for x in indices]
#overlap_red = [overlap[x] for x in indices]
