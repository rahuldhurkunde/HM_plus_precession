import h5py as hf
import numpy as np
import matplotlib.pyplot as plt
import functions as func
import sys

def plot_tau_function(signal_tau, tau_diff):
	tau_tolerance = 0.1
	bin_edges, statistic = func.compute_tauThreshold_envelope(signal_tau, tau_diff, 50)  #50 is number of bins for the histogram binning
	tau_func = func.fit_tau_envelope(bin_edges, statistic, tau_tolerance)
	#plt.plot(signal_tau, tau_diff, 'x')
	#plt.show()
	return bin_edges, statistic

nsplits = 200 

dir_name = sys.argv[1] 

tau_diff = np.array([])
signal_tau = np.array([]) 
for k in range(nsplits):
    filename = "%s/%s/output/TF2-FF-0-10.hdf" %(dir_name, k)
    f = hf.File(filename, 'r')
    tau_diff = np.concatenate([tau_diff, abs(np.array(f['tau0']) - np.array(f['rec_tau']))])
    signal_tau = np.concatenate([signal_tau, np.array(f['tau0'])])

bin_edges, statistic = plot_tau_function(signal_tau, tau_diff)

filename = "/work/rahul.dhurkunde/HM_and_precession/injections/tau_files/tau_crawl_50000_%s.txt" %dir_name
np.savetxt(filename, np.c_[bin_edges, statistic])
