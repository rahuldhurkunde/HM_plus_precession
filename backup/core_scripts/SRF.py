#python SRF.py --compute_SRF 1 --inj_class 0 --inj_dir /work/rahul.dhurkunde/HM_and_precession/injections/100000_inj/aligned_injections --approximant IMRPhenomD --HMs 0 --psd_file /work/rahul.dhurkunde/HM_and_precession/psds/ZDHP.txt --nsets 200

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

def read_injections(nsets, nsignals, f_min, inj_dir):
	sg = []
	for k in range(nsets):
		filename = '%s/%s.hdf' %(inj_dir, k)
		#print('loaded', filename)
		sg = sg + inj.read_injections_HDF(filename, nsignals, f_min)

	sg_m1 = [x.m1 for x in sg]
	sg_m2 = [x.m2 for x in sg]
	mtotal = sg_m1 + sg_m2
	q = np.divide(sg_m1, sg_m2)
	end = time.time()
	return sg, sg_m1, sg_m2, mtotal, q

def read_FFs(inj_type, nsets):
	FF = np.array([])
	m1 = np.array([])
	m2 = np.array([])
	for k in range(nsets):
		filename = "../%s/output/H1-FF_%s-0-10.hdf" %(inj_type, k)
		f = h5py.File(filename, 'r')

		FF = np.concatenate([FF, np.array(f['FF'])])
		m1 = np.concatenate([m1, np.array(f['mass1'])])
		m2 = np.concatenate([m2, np.array(f['mass2'])])
	return FF

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

def save_SRF(savefile, xedge, yedge, SRF):
	with h5py.File(savefile, 'w') as f:
		f.create_dataset("SRF", data=SRF)
		f.create_dataset("xedge", data=xedge)
		f.create_dataset("yedge", data=yedge)

def save_injections_inside(nbinsx, nbinsy, xedge, yedge, indices, sigmasq_list):
	for y in range(nbinsy):
		for x in range(nbinsx):
			binno = x + y*nbinsx
			savefile = "SRF_results/%sx%s/indices/%s.hdf" %(nbinsx, nbinsy, binno)
			x1 = np.full(20, xedge[x])
			x2 = np.full(20, xedge[x+1])
			y1 = np.full(20, yedge[y])
			y2 = np.full(20, yedge[y+1])
			with h5py.File(savefile, 'w') as f:
				f.create_dataset("ind",data=indices[binno])
				f.create_dataset("x1", data=x1)
				f.create_dataset("x2", data=x2)
				f.create_dataset("y1", data=y1)
				f.create_dataset("y2", data=y2)
				f.create_dataset("sigmasq", data=sigmasq_list[binno])
				

def get_vmax_vmin(prev_vmax, prev_vmin, SRF):
	if (max(SRF) > prev_vmax):
		vmax = max(SRF)
	elif (min(SRF) < prev_vmin):
		vmin = min(SRF) 
	else:
		vmax = prev_vmax
		vmin = prev_vmin
	return vmax, vmin

def create_bins(sg_m1, sg_m2, nbinsx, nbinsy):
	xbins = np.linspace(min(sg_m1)+0.5, max(sg_m1)-0.5, nbinsx+1)
	ybins = np.linspace(min(sg_m2)+0.1, max(sg_m2)-0.1, nbinsy+1)
	#xvalues, yvalues = np.meshgrid(xbins, ybins)
	return xbins, ybins

def find_injIndices_inside_bins(xparameter, yparameter, xbins, ybins):
	injIndices = []
	xcells = np.digitize(xparameter, xbins)
	ycells = np.digitize(yparameter, ybins)
	xcells = xcells - 1
	ycells = ycells - 1
	for k in range(len(sg)):
		cell_index = int(ycells[k]*(len(xbins)-1) + xcells[k])
		injIndices.append(np.array([xparameter[k], yparameter[k], xcells[k], ycells[k], cell_index]))
	return injIndices

def get_corners_from_cell_ind(xbins, ybins, cell_no):
	col = int(cell_no%(len(xbins)-1))
	row = int((cell_no - col)/(len(xbins)-1))
	x1 = col 
	x2 = col + 1
	y1 = row
	y2 = row + 1
	return x1, y1, x2, y2 

def check_bin_cell_indices(injIndices, xbins, ybins, cell_no):
	x1, y1, x2, y2 = get_corners_from_cell_ind(xbins, ybins, cell_no)
	plt.plot(xbins[x1], ybins[y1], 'x', color='r')
	plt.plot(xbins[x1], ybins[y2], 'x', color='r')
	plt.plot(xbins[x2], ybins[y1], 'x', color='r')
	plt.plot(xbins[x2], ybins[y2], 'x', color='r')
	plt.xlim([min(sg_m1), max(sg_m1)])
	plt.ylim([min(sg_m2), max(sg_m2)])

	in_indices = get_injs_inside_cell(injIndices, cell_no)
	print('Total injections inside', len(in_indices), "Corner indices are", xbins[x1], ybins[y1], xbins[x2], ybins[y2])
	for k in range(len(in_indices)):
		plt.scatter(sg_m1[in_indices[k]], sg_m2[in_indices[k]])
	plt.show()	

def get_injInside(binnumber, bin_x, bin_y):
	inj_with_row = np.where(binnumber[0] == bin_x+1) 
	inj_with_col = np.where(binnumber[1] == bin_y+1) 
	injInd = np.intersect1d(inj_with_row, inj_with_col)
	return injInd

def compute_SRF(param_x, param_y, nbinsx, nbinsy, FF, sg, approximant, psd, f_min, delta_f, detector, HMs):
	start = time.time()
	ret = binned_statistic_2d(param_x, param_y, FF, statistic='count', bins=[nbinsx, nbinsy], expand_binnumbers=True)
	end = time.time()
	print('Binned statistics computed in', end-start)

	SRF = []
	indices = []
	sigmasq_list = []
	ncells = int(nbinsx*nbinsy)
	binnumber = ret.binnumber
	count = ret.statistic
	xedge = ret.x_edge
	yedge = ret.y_edge
	print('MKC', xedge[0], yedge[0], min(param_x), min(param_y)) 

	for y in range(nbinsy):
		for x in range(nbinsx):
			numerator = 0.0
			denominator = 0.0
			injInd = np.array([])
			ninjs = int(count[x, y])
			sigma = []

			if(ninjs !=0):	
				injInd = get_injInside(binnumber, x, y)
				for n in range(ninjs):
					inj_ind = injInd[n]
					if HMs==True:
						sp, sc = func.generate_signal(sg[inj_ind], delta_f, f_min, approximant)
					else:
						modes = [[2,2],[2,-2]]
						sp, sc = func.generate_signal(sg[inj_ind], delta_f, f_min, approximant, modes=modes)
					s_f = func.signal_to_detector_frame(detector, sp, sc, sg[inj_ind])
					s_f.resize(len(psd))
					temp_overlap = pycbc.filter.matchedfilter.overlap(s_f, s_f, psd=psd, low_frequency_cutoff=f_min, normalized=False)
					numerator += FF[inj_ind]**3*temp_overlap**3
					denominator += temp_overlap**3
					sigma.append(temp_overlap)
				temp_SRF = numerator/denominator
			else:
				temp_SRF = np.nan
				sigma.append(0.0)
				print('No injections found in bin', x, y)

			print([x,y], 'Inj inside', len(injInd), 'SRF', temp_SRF)
			indices.append(injInd)
			SRF.append(temp_SRF)
			sigmasq_list.append(sigma)
	return SRF, indices, sigmasq_list, xedge, yedge

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
	#im = ax.scatter(xx, yy, c = SRF, marker='s', vmax=vmax, vmin=vmin, cmap = 'YlGnBu')
	#im = ax.scatter(mtotal, q, c = SRF, marker='s', vmax=vmax, vmin=vmin, cmap = 'YlGnBu')
	return im

def combined_plot(SRF_array, inj_types, xbins, ybins, fig, axs, vmax, vmin):
	for row in range(2):
		for col in range(2):
			current_inj = inj_types[row*2 + col]
			print(col+row*2, current_inj)
			ax = axs[row, col]
			#ax = axs
			
			FF = read_FFs(current_inj, nsets)
			SRF = SRF_array[col + row*2]

			vmax, vmin = get_vmax_vmin(vmax, vmin, SRF)
			print(current_inj, 'After smoothing SRF range', min(SRF), max(SRF))
			
			im = plot_SRF(SRF, xbins, ybins, fig, ax, vmax, vmin)
			ax.set_title(current_inj)

	fig.subplots_adjust(right=0.8, hspace = 0.3)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
	fig.colorbar(im, cax=cbar_ax)
	fig.text(0.5, 0.04, '$m_1$', ha='center')
	fig.text(0.04, 0.5, '$m_2$', va='center', rotation='vertical')
	fig.text(0.5, 0.95, 'SRF for ZDHP', ha='center', fontsize='x-large')	
	#plt.savefig("zdhp_SRF.png", dpi=600)
	plt.show()
	

f_min = 30.0
fref = 100.0
pref = 0.0
delta_f = 1.0/32
delta_t = 1.0/2048
detector = detector.Detector('H1')
nsignals = 10
nbinsx = 12
nbinsy = 10 
fig, axs = plt.subplots(2,2)
inj_types = ['aligned', 'onlyHM', 'onlyPrecession', 'HM_prec']

#Command line parser
parser = ap.ArgumentParser()
parser.add_argument('--compute_SRF',
			type=int, default = 0,
			help = 'Specify 1 to compute SRF (default option will only plot the SRF)')
parser.add_argument('--inj_class',
			type = int,
			help = 'Injection class (integer) 0-aligned, 1-onlyHM, 2-onlyPrecession, 3-HM_prec')
parser.add_argument('--inj_dir',
			help = 'Specify injections dir')
parser.add_argument('--psd_file',
			help = 'Specify the ASD file')
parser.add_argument('--approximant',
			help = 'Specify approximant accordingly to the Injection class')
parser.add_argument('--HMs',
			type=int,
			help = 'Specify HMs True=1 or False=0 accordingly to the Injection class')
parser.add_argument('--nsets',
			type = int, required=True,
			help = 'Specify number of injection files')
	
args, remaining_args = parser.parse_known_args()
nsets = args.nsets 

if args.compute_SRF == 1:
	print("Computing SRF", args.compute_SRF)
	srf_dir = 'SRF_results/%sx%s' %(nbinsx, nbinsy)
	indices_dir = 'SRF_results/%sx%s/indices' %(nbinsx, nbinsy)
	if(os.path.isdir(srf_dir) == False):
		raise ValueError('Create a SRF dir')
	if(os.path.isdir(indices_dir) == False):
		raise ValueError('Create a /indices dir inside SRF dir')

	current_inj = inj_types[args.inj_class]
	if(args.HMs == 0):
		HMs = False
	else:
		HMs = True
	print(current_inj, "=", args.approximant, HMs, args.inj_dir, 'PSD =', args.psd_file, 'nsets', nsets)

	length = int(1.0/(delta_f*delta_t)/2 + 1)
	psd = pycbc.psd.read.from_txt(args.psd_file, length, delta_f, f_min, is_asd_file=True)
	
	start = time.time()
	sg, sg_m1, sg_m2, mtotal, q = read_injections(nsets, nsignals, f_min, args.inj_dir)
	end = time.time()
	print("Inj loading time", end-start)

	start = time.time()
	FF = read_FFs(current_inj, nsets)
	SRF, indices, sigmasq_list, xedge, yedge = compute_SRF(sg_m1, sg_m2, nbinsx, nbinsy, FF, sg, args.approximant, psd, f_min, delta_f, detector, HMs)
	end = time.time()
	print("Time taken for SRF computation", end-start)

	SRFfile = 'SRF_results/%sx%s/%s_%s.hdf' %(nbinsx, nbinsy, current_inj, nsets)
	print('Saving SRF in file', SRFfile)
	save_SRF(SRFfile, xedge, yedge, SRF)
	save_injections_inside(nbinsx, nbinsy, xedge, yedge, indices, sigmasq_list)

else:
	#Plot the SRF
	xbins, ybins, SRF = read_SRF(nbinsx, nbinsy, inj_types, nsets)
	vmax = np.nanmax(SRF)
	vmin = np.nanmin(SRF)
	#vmin = 0.1
	print('Vmax', vmax, 'Vmin', vmin)

	combined_plot(SRF, inj_types, xbins, ybins, fig, axs, vmax, vmin)

