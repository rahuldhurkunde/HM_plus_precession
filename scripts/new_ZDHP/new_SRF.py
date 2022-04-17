import numpy as np
import matplotlib.pyplot as plt
import h5py 
import time
import injections as inj
import functions as func
import pycbc
from scipy.interpolate import griddata
from scipy.stats import binned_statistic_2d
from pycbc import psd, detector
import lal
import seaborn as sns
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import pandas as pd
import argparse as ap
from matplotlib.patches import Rectangle
import os
lal.ClobberDebugLevel(0)

def read_FFs(inj_type):
	filename = "%s/10000_FFs.hdf" %(inj_type)
	f = h5py.File(filename, 'r')

	FF = np.array(f['FF'])
	m1 = np.array(f['mass1'])
	m2 = np.array(f['mass2'])
	sigmasqs = np.array(f['sigmasqs'])
	return FF, sigmasqs, m1, m2

def read_SRF(nbinsx, nbinsy, inj_type, nsets):
	SRF = [] 
	for row in range(2):
		for col in range(2):
			current_inj = inj_type[row*2 + col]
			srf_file = "SRF_results/%sx%s/%s_%s.hdf" %(nbinsx, nbinsy, current_inj, nsets)
			hf = h5py.File(srf_file)
			SRF.append(hf['SRF'])
		xedge = hf['xedge']
		yedge = hf['yedge']
	return xedge, yedge, np.array(SRF)

def get_injInside(binnumber, bin_x, bin_y):
	inj_with_row = np.where(binnumber[0] == bin_x+1) 
	inj_with_col = np.where(binnumber[1] == bin_y+1) 
	injInd = np.intersect1d(inj_with_row, inj_with_col)
	return injInd

def compute_SRF(param_x, param_y, nbinsx, nbinsy, FF, sigmasqs):
	start = time.time()
	ret = binned_statistic_2d(param_x, param_y, FF, statistic='count', bins=[nbinsx, nbinsy], expand_binnumbers=True)
	end = time.time()
	print('Binned statistics computed in', end-start)

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

			print([x,y], 'Inj inside', len(injInd), 'SRF', temp_SRF)
			indices.append(injInd)
			SRF.append(temp_SRF)
	return np.array(SRF), indices, xedge, yedge

def plot_SRF(SRF, xedge, yedge, fig, ax, vmax, vmin):
	levels = MaxNLocator(nbins=15).tick_values(vmin, vmax)
	cmap = plt.cm.YlGnBu
	norm = BoundaryNorm(levels, ncolors = cmap.N, clip=True)	

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
	im = ax.pcolormesh(xx, yy, z, cmap='YlGnBu', norm=norm, shading='auto')	
	return im

nbinsx = 10 
nbinsy = 8 
fig, axs = plt.subplots(1,1)
inj_types = ['aligned', 'onlyHM', 'onlyPrecession', 'HM_prec']
fig, ax = plt.subplots(1,1)

#Command line parser
parser = ap.ArgumentParser()
parser.add_argument('--compute_SRF',
			type=int, default = 0,
			help = 'Specify 1 to compute SRF (default option will only plot the SRF)')
parser.add_argument('--inj_class',
			type = int, required=True,
			help = 'Injection class (integer) 0-aligned, 1-onlyHM, 2-onlyPrecession, 3-HM_prec')
	
args, remaining_args = parser.parse_known_args()

if args.compute_SRF == 1:
	print("Computing SRF for bins", nbinsx, 'x', nbinsy)
	current_inj = inj_types[args.inj_class]
	print(current_inj)
	
	start = time.time()
	FF, sigmasqs, sg_m1, sg_m2 = read_FFs(current_inj)
	SRF, indices, xedge, yedge = compute_SRF(sg_m1, sg_m2, nbinsx, nbinsy, FF, sigmasqs)
	end = time.time()
	print("Time taken for SRF computation", end-start)

	vmax = np.nanmax(SRF)
	vmin = np.nanmin(SRF)
	#vmin = 0.1
	print('Vmax', vmax, 'Vmin', vmin)

	plot_SRF(SRF, xedge, yedge, fig, ax, vmax, vmin)
	plt.show()
