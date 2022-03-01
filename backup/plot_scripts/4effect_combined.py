import h5py as hf
import numpy as np
import matplotlib.pyplot as plt
import scipy
import injections as inj
import functions as func
from scipy.interpolate import griddata
import scipy.ndimage as ndimage
from scipy.stats import binned_statistic_2d
import lalsimulation as lalsim
import pycbc
from pycbc import waveform, psd, filter, conversions
import lal
import time
lal.ClobberDebugLevel(0)

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

def read_injections(inj_dir, nsets, nsignals, f_min, fref, pref):
	sg = []
	for k in range(nsets):
		filename = '%s/%s.hdf' %(inj_dir, k)
		#print('loaded', filename)
		sg = sg + inj.read_injections_HDF(filename, nsignals, f_min)

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
	chi_p = np.array(conversions.chi_p(sg_m1[:], sg_m2[:], spin1x[:], spin1y[:], spin2x[:], spin2y[:]))

	thetaJN = get_thetaJN(sg_m1, sg_m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, inclination, fref, pref)
	inclination = np.degrees(inclination)
	return sg, sg_m1, sg_m2, chi_p, thetaJN, inclination

def binned_plot(x, y, z, xbins, ybins, fig, ax, vmin, vmax):
	x_bins = np.linspace(min(x), max(x), xbins)
	y_bins = np.linspace(min(y), max(y), ybins)

	ret = binned_statistic_2d(x, y, z, statistic = np.mean, bins=[x_bins, y_bins])
	print('After smoothing max value', np.nanmax(ret.statistic), 'min value', np.nanmin(ret.statistic))
	im = ax.imshow(ret.statistic.T, origin='lower', extent=[x_bins[0], x_bins[-1], y_bins[0], y_bins[-1]], vmin=vmin, vmax=vmax, aspect='auto', cmap = 'YlGnBu')
	#im = ax.imshow(ret.statistic.T, origin='lower', vmin=vmin, vmax=vmax)
	#ax.set_xlim([min(x), max(x)])
	#ax.set_ylim([min(y), max(y)])
	return im

def nonuniform_imshow(x, y, z, npoints, fig, ax, levels):
	# Create regular grid
	xi, yi = np.linspace(x.min(), x.max(), npoints), np.linspace(y.min(), y.max(), npoints)
	xi, yi = np.meshgrid(xi, yi)

  	# Interpolate missing data
	zi = griddata((x, y), z, (xi, yi), method='linear')
	Z2 = ndimage.gaussian_filter(zi, sigma=0.5, order=0)
	im = ax.contourf(xi, yi, Z2, 15, levels=levels)
	return im 

def plot_cdf(FF,fig_titles):
	X = np.sort(FF)
	Y = np.array(range(len(FF)))/float(len(FF))
	plt.yscale('log')
	plt.plot(X,Y, label=fig_titles)
	#plt.grid()
	plt.xlabel('Fitting factors')
	plt.ylabel('Cummulative distribution')
	#plt.show()

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
directories = ['aligned', 'onlyHM', 'onlyPrecession', 'HM_prec']
#directories = [ 'onlyHM', 'onlyPrecession']
fig_titles = ['aligned_spins_noHM', 'aligned_HM', 'Precessing_noHM', 'HM_and_precession']
injection_dir = ['/work/rahul.dhurkunde/HM_and_precession/injections/50000_inj/aligned_injections', '/work/rahul.dhurkunde/HM_and_precession/injections/50000_inj/nonaligned_injections']

vmax = 1.0                  
vmin = 0.76
fref = 100.0
pref = 0.0
f_min = 30.0
levels = np.linspace(vmin, vmax, 8)
nsignals = 100
nsets = 1000
xbins = 20
ybins = 20

fig, axs = plt.subplots(1, 2, figsize=(12, 22))
for row in range(1):
	for col in range(2):
		#inj_dir = injection_dir[col]
		#print('Loading injections from', inj_dir)
		#start = time.time()
		#sg, sg_m1, sg_m2, chi_p, thetaJN, inclination = read_injections(inj_dir, nsets, nsignals, f_min, fref, pref)
		#end = time.time()
		#print("Inj loading time", end-start, 'Total injections', len(sg))

		current_dir = directories[row*2 + col]
		print(col+row, current_dir) 
		#ax = axs[row, col]		
		ax = axs[col]

		FF = np.array([])
		tau_diff = np.array([])
		m1 = np.array([])
		m2 = np.array([])
		for k in range(nsets):
			filename = "../%s/%s/output/TF2-FF-0-10.hdf" %(current_dir, k)
			f = hf.File(filename, 'r')

			FF = np.concatenate([FF, np.array(f['FF'])])
			tau_diff = np.concatenate([tau_diff, abs(np.array(f['tau0']) - np.array(f['rec_tau']))])
			m1 = np.concatenate([m1, np.array(f['mass1'])])
			m2 = np.concatenate([m2, np.array(f['mass2'])])

		mtotal = m1 + m2
		q = np.divide(m1, m2)
		#im = ax.scatter(m1, m2, c = FF, vmax=1.0, vmin=0.6)
		#im = nonuniform_imshow(m1, m2, FF, 100, fig, ax, levels)
		im = binned_plot(m1, m2, FF, xbins, ybins, fig, ax, vmin, vmax)
		#im = binned_plot(q, inclination, FF, xbins, ybins, fig, ax, vmin, vmax)
		ax.set_title(fig_titles[col+row*2])
		#plot_cdf(FF, fig_titles[col+row*2])

print("Min FF", min(FF))
fig.subplots_adjust(right=0.8, hspace = 0.3)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)
fig.text(0.5, 0.04, '$q$', ha='center')
fig.text(0.04, 0.5, 'inclination', va='center', rotation='vertical')
fig.text(0.5, 0.95, 'FF for ZDHP', ha='center', fontsize='x-large')
#plt.savefig('smooth_4FF_m1_m2.png')

#plt.legend()
#plt.grid()
#plt.savefig('4cdf.png', dpi = 600)
plt.show()
