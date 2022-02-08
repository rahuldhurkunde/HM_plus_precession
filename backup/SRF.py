import numpy as np
import matplotlib.pyplot as plt
import h5py 
import time
import injections as inj
import functions as func
import pycbc
from pycbc import psd, detector
import lal
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
		filename = "../%s/%s/output/TF2-FF-0-10.hdf" %(inj_type, k)
		f = h5py.File(filename, 'r')

		FF = np.concatenate([FF, np.array(f['FF'])])
		m1 = np.concatenate([m1, np.array(f['mass1'])])
		m2 = np.concatenate([m2, np.array(f['mass2'])])
	return FF

def read_SRF(inj_type, nsets):
	SRF = [] 
	for row in range(2):
		if row == 0:
			col_range = 1
		else:
			col_range = 2
		for col in range(2):
			current_inj = inj_type[row*2 + col]
			srf_file = "zdhp_%s_%s.txt" %(current_inj, nsets)
			temp = np.array(np.loadtxt(srf_file))
			SRF.append(temp)
	return np.array(SRF)

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

def find_injIndices_inside_bins(sg, xbins, ybins):
	injIndices = []
	xcells = np.digitize(sg_m1, xbins)
	ycells = np.digitize(sg_m2, ybins)
	xcells = xcells - 1
	ycells = ycells - 1
	for k in range(len(sg)):
		cell_index = int(ycells[k]*(len(xbins)-1) + xcells[k])
		injIndices.append(np.array([sg_m1[k], sg_m2[k], xcells[k], ycells[k], cell_index]))
	return injIndices

def get_corners_from_cell_ind(xbins, ybins, cell_no):
	col = int(cell_no%(len(xbins)-1))
	row = int((cell_no - col)/(len(xbins)-1))
	x1 = col 
	x2 = col + 1
	y1 = row
	y2 = row + 1
	return x1, y1, x2, y2 

def get_injs_inside_cell(injIndices, cell_index):     #Search can be done in a Pythonic way
	in_indices = []
	for k in range(len(injIndices)):
		if(injIndices[k][4]==cell_index):
			in_indices.append(k)
	return(np.array(in_indices))

def check_bin_cell_indices(injIndices, xbins, ybins, cell_no):
	x1, y1, x2, y2 = get_corners_from_cell_ind(xbins, ybins, cell_no)
	print(cell_no, "Corner indices are", x1, y1, x2, y2)
	plt.plot(xbins[x1], ybins[y1], 'x', color='r')
	plt.plot(xbins[x1], ybins[y2], 'x', color='r')
	plt.plot(xbins[x2], ybins[y1], 'x', color='r')
	plt.plot(xbins[x2], ybins[y2], 'x', color='r')
	plt.xlim([min(sg_m1), max(sg_m1)])
	plt.ylim([min(sg_m2), max(sg_m2)])

	in_indices = get_injs_inside_cell(injIndices, cell_no)
	for k in range(len(in_indices)):
		plt.scatter(sg_m1[in_indices[k]], sg_m2[in_indices[k]])
	plt.show()	

def compute_SRF(xbins, ybins, FF, sg, injIndices, approximant, psd, f_min, delta_f, detector, HMs):
	ncells = (len(xbins)-1)*(len(ybins)-1) 
	print("Total cells", ncells)
	SRF = []
	for k in range(ncells):
		in_indices = get_injs_inside_cell(injIndices, k)
		if (len(in_indices) != 0):
			numerator = 0.0
			denominator = 0.0
			for n in range(len(in_indices)):
				inj_ind = in_indices[n]
				if HMs==True:
					sp, sc = func.generate_signal(sg[inj_ind], delta_f, f_min, approximant)
				else:
					modes = [[2,2],[2,-2]]
					sp, sc = func.generate_signal(sg[inj_ind], delta_f, f_min, approximant, modes=modes)
				s_f = func.signal_to_detector_frame(detector, sp, sc, sg[inj_ind])
				s_f.resize(len(psd))
				overlap = pycbc.filter.matchedfilter.overlap(s_f, s_f, psd=psd, low_frequency_cutoff=f_min, normalized=False)
				numerator += FF[inj_ind]**3*overlap**3
				denominator += overlap**3
			temp_SRF = numerator/denominator	
		else:
			temp_SRF = np.nan 
		print(k, "Cell number done") 
		SRF.append(temp_SRF)
	return SRF

def plot_SRF(SRF, xbins, ybins, fig, ax, vmax, vmin):
	temp_x = xbins[:-1]
	temp_y = ybins[:-1]
	delta_x = xbins[1]-xbins[0]
	delta_y = ybins[1]-ybins[0]
	x_new = temp_x + delta_x/2.0
	y_new = temp_y + delta_y/2.0
	xx, yy = np.meshgrid(x_new, y_new)
	im = ax.scatter(xx, yy, c = SRF, vmax=vmax, vmin=vmin)
	return im

def combined_plot(SRF, inj_types, xbins, ybins, fig, axs, vmax, vmin):
	for row in range(2):
		if row == 0:
			col_range = 1
		else:
			col_range = 2
		for col in range(2):
			current_inj = inj_types[row*2 + col]
			print(col+row, current_inj)
			ax = axs[row, col]
			
			FF = read_FFs(current_inj, nsets)
			srf_file = "zdhp_%s_%s.txt" %(current_inj, nsets)
			SRF = np.loadtxt(srf_file)

			vmax, vmin = get_vmax_vmin(vmax, vmin, SRF)
			print(current_inj, 'After smoothing SRF range', min(SRF), max(SRF))
			
			im = plot_SRF(SRF, xbins, ybins, fig, ax, vmax, vmin)
			ax.set_title(current_inj)

	fig.subplots_adjust(right=0.8, hspace = 0.3)
	cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
	fig.colorbar(im, cax=cbar_ax)
	fig.text(0.5, 0.04, '$m_1$', ha='center')
	fig.text(0.04, 0.5, '$m_2$', va='center', rotation='vertical')
	plt.show()
	

fig, axs = plt.subplots(2,2)
inj_types = ['aligned', 'onlyHM', 'onlyPrecession', 'HM_prec']

f_min = 30.0
fref = 100.0
pref = 0.0
delta_f = 1.0/32
delta_t = 1.0/4096
detector = detector.Detector('H1')
approximant = 'IMRPhenomXPHM'

HMs = True
inj_dir = '/work/rahul.dhurkunde/HM_and_precession/injections/50000_inj/nonaligned_injections'
psd_file = '/work/rahul.dhurkunde/HM_and_precession/psds/ZDHP.txt'
nsets = 200 
nsignals = 100
nbinsx = 25
nbinsy = 15
vmax = 1.0
vmin = 0.3

length = int(1.0/(delta_f*delta_t)/2 + 1)
psd = pycbc.psd.read.from_txt(psd_file, length, delta_f, f_min, is_asd_file=True)

start = time.time()
sg, sg_m1, sg_m2, mtotal, q = read_injections(nsets, nsignals, f_min, inj_dir)
end = time.time()
print("Inj loading time", end-start)
xbins, ybins = create_bins(sg_m1, sg_m2, nbinsx, nbinsy)
injIndices = find_injIndices_inside_bins(sg, xbins, ybins)

SRF = read_SRF(inj_types, nsets)
vmax = np.amax(SRF) 
vmin = np.amin(SRF)

combined_plot(SRF, inj_types, xbins, ybins, fig, axs, vmax, vmin)

start = time.time()
#SRF = compute_SRF(xbins, ybins, FF, sg, injIndices, approximant, psd, f_min, delta_f, detector, HMs)
#np.savetxt('zdhp_both_200.txt', SRF)
end = time.time()
print("Time taken for SRF computation", end-start)

#for i in range(8):
#	check_bin_cell_indices(injIndices, xbins, ybins, i)
