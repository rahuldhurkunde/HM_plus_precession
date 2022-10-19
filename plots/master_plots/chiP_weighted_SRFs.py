import numpy as np
import matplotlib.pyplot as plt
import h5py 
import time
import injections as inj
import functions as func
import pycbc
from scipy.stats import binned_statistic_2d
from pycbc import psd, detector, conversions
import lal
import lalsimulation as lalsim
from matplotlib import colors
from matplotlib.colors import LinearSegmentedColormap
import os
lal.ClobberDebugLevel(0)
plt.rcParams.update({
    "text.usetex": True})
from matplotlib.ticker import AutoMinorLocator, FormatStrFormatter
lalParams = lal.CreateDict()

my_gradient = LinearSegmentedColormap.from_list('my_gradient', (
    # Edit this gradient at https://eltos.github.io/gradient/#9.8:FB7D28-31.2:F2DE92-67.5:F6F6F6-77.6:B8DEE0-100:3456DC
    (0.000, (0.984, 0.490, 0.157)),
    (0.098, (0.984, 0.490, 0.157)),
    (0.312, (0.949, 0.871, 0.573)),
    (0.675, (0.965, 0.965, 0.965)),
    (0.776, (0.722, 0.871, 0.878)),
    (1.000, (0.204, 0.337, 0.863))))

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

def read_injections(inj_dir, fref, pref):
	filename = "/work/rahul.dhurkunde/HM_and_precession/injections/final_injections/%s/combined_inj.hdf" %(inj_dir)
	hf = h5py.File(filename, 'r')
	
	sg_m1 = hf['mass1'][:]
	sg_m2 = hf['mass2'][:]
	spin1x = hf['spin1x'][:]
	spin1y = hf['spin1y'][:]
	spin1z = hf['spin1z'][:]
	spin2x = hf['spin2x'][:]
	spin2y = hf['spin2y'][:]
	spin2z = hf['spin2z'][:]
	inclination = hf['inclination'][:]

	chi_p = np.array(conversions.chi_p(sg_m1[:], sg_m2[:], spin1x[:], spin1y[:], spin2x[:], spin2y[:]))
	thetaJN = get_thetaJN(sg_m1, sg_m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, inclination, fref, pref)
	inclination = np.degrees(inclination)
	return sg_m1, sg_m2, chi_p, thetaJN, inclination

def read_FFs(filename):
	f = h5py.File(filename, 'r')

	FF = np.array(f['FF'])
	m1 = np.array(f['mass1'])
	m2 = np.array(f['mass2'])
	sigmasqs = np.array(f['sigmasqs'])
	return FF, sigmasqs, m1, m2

def get_injInside(binnumber, bin_x, bin_y):
    inj_with_row = np.where(binnumber[0] == bin_x+1)
    inj_with_col = np.where(binnumber[1] == bin_y+1)
    injInd = np.intersect1d(inj_with_row, inj_with_col)
    return injInd

def compute_SRF(param_x, param_y, nbinsx, nbinsy, FF, sigmasqs):
	start = time.time()
	ret = binned_statistic_2d(param_x, param_y, FF, statistic='count', bins=[nbinsx, nbinsy], expand_binnumbers=True)
	end = time.time()

	SRF = []
	indices = []
	ncells = int(nbinsx*nbinsy)
	binnumber = ret.binnumber
	count = ret.statistic
	xedge = ret.x_edge
	yedge = ret.y_edge

	for y in range(nbinsy):
		for x in range(nbinsx):
			numerator = 0.0
			denominator = 0.0
			injInd = np.array([])
			ninjs = int(count[x, y])

			if(ninjs !=0):	
				injInd = get_injInside(binnumber, x, y)
				for n in range(ninjs):
					inj_ind = injInd[n]
					numerator += FF[inj_ind]**3*sigmasqs[inj_ind]**3
					denominator += sigmasqs[inj_ind]**3
				temp_SRF = numerator/denominator
			else:
				temp_SRF = np.nan
				print('No injections found in bin', x, y)

			#print([x,y], 'Inj inside', len(injInd), 'SRF', temp_SRF)
			indices.append(injInd)
			SRF.append(temp_SRF)
	return np.array(SRF), indices, xedge, yedge

def plot_SRF(SRF, xedge, yedge, fig, ax, vmax, vmin):
	#levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
	cmap = plt.cm.YlGnBu
	#norm = BoundaryNorm(levels, ncolors = cmap.N, clip=True)	

	temp_x = xedge[:-1]
	temp_y = yedge[:-1]
	dx = xedge[1]-xedge[0]
	dy = yedge[1]-yedge[0]
	x_new = temp_x + dx/2.0
	y_new = temp_y + dy/2.0
	xx, yy = np.meshgrid(x_new, y_new)

	mtotal = xx + yy
	q = np.divide(xx, yy)
	z = SRF.reshape(len(y_new), len(x_new))
	#z[0,-1] = np.nan 
	im = ax.pcolormesh(xx, yy, z, cmap=my_gradient, vmax=vmax, vmin=vmin, shading='nearest') 
	return im


ref_SNR = [9.37, 9.92]
thresh_SNR = [9.37, 9.70, 10.44, 10.44]
fref = 100.0
pref = 0.0

nbinsx = 8 
nbinsy = 7 
fig, axs = plt.subplots(2,2)

psd_list = ['ZDHP', 'aplus', 'voyager']
fig_title = ['aLIGO design', 'A+', 'Voyager']
	
sg_m1, sg_m2, chi_p, thetaJN, inclination =	read_injections('nonaligned_injections', fref, pref)
mtotal = sg_m1 + sg_m2
q = np.divide(sg_m1, sg_m2)

xparam = q
yparam = chi_p

for row in range(2):
	for col in range(2):
		if (row*2 + col >= 3):
			break
		current_psd = psd_list[col + row*2]
		filename = '../%s/HM_prec/50000_FFs.hdf' %current_psd
		hf = h5py.File(filename, 'r')
		ax = axs[row, col]
	
		FF, sigmasqs, sg_m1, sg_m2 = read_FFs(filename)
		SRF, indices, xedge, yedge = compute_SRF(xparam, yparam, nbinsx, nbinsy, FF, sigmasqs)

		weighted_SRF = SRF*(thresh_SNR[2]/ref_SNR[1])**3

		vmax = np.nanmax(weighted_SRF)
		vmin = np.nanmin(weighted_SRF)
		print('Vmax', vmax, 'Vmin', vmin)
		im = plot_SRF(weighted_SRF, xedge, yedge, fig, ax, 1.13, 0.74)
		ax.set_title(fig_title[col + row*2])

fig.subplots_adjust(right=0.8, hspace = 0.3)
cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
fig.colorbar(im, cax=cbar_ax)

#plt.savefig('background_adjusted_SRF_diffPSDs.jpg', dpi=600)
plt.show()	
