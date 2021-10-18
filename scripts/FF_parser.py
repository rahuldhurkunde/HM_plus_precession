import argparse as ap
import configparser
import argcomplete
import functions as func
import injections as inj
import pycbc
from pycbc import waveform, conversions, filter, types, distributions, detector, psd
import matplotlib.pyplot as plt
import numpy as np
import mass as mass
import functions as func
import injections as inj
import scipy
import seaborn as sns
import importlib
from random import randrange
import time
import matplotlib.gridspec as gridspec

def determine_values_from_parsing_config(config):
	f_min = config.getfloat('Required', 'f_min')
	delta_f = 1.0/(config.getfloat('Required', 'sampling_freq'))
	delta_t = 1.0/(config.getfloat('Required', 'sampling_rate'))
	nsignal = config.getint('Required', 'nsignal')
	psdname = config['Required'][ 'psd']
	detectorname = config['Required']['detector']
	nbins_tau = config.getfloat('Tau_threshold', 'nbins_tau')
	tau_tolerance = config.getfloat('Tau_threshold', 'tau_tolerance')
	return f_min, delta_f, delta_t, psdname, detectorname, nsignal, nbins_tau, tau_tolerance

parser = ap.ArgumentParser()
parser.add_argument('-c', '--config', required=True, help='specify config file path')
parser.add_argument('-tb', '--template_bank', help='specify template bank', default = "../banks/vanilla_BBHbank")
parser.add_argument('-inj', '--injections', help='specify injection file', default = "../injections/random_injections.hdf5")
parser.add_argument('--tau_threshold', help='specify tau thresholds file', default = "../injections/tau_threshold.txt")

args, remaining_argv = parser.parse_known_args()

if args.config:
	config = configparser.ConfigParser()
	config.read([args.config])

f_min, delta_f, delta_t, psdname, detectorname, nsignal, nbins_tau, tau_tolerance = determine_values_from_parsing_config(config)

#Read template bank
tb = func.read_tb(args.template_bank, f_min)

#Read injections
sg = inj.read_injections_HDF(args.injections, nsignal)

#Tau_threshold function
tau_bin_edges = np.loadtxt(args.tau_threshold)[:,0]
tau_bin_statistic = np.loadtxt(args.tau_threshold)[:,1]
tau_func = func.fit_tau_envelope(tau_bin_edges, tau_bin_statistic, tau_tolerance)

#Read/compute matches
length = int(1.0/(delta_f*delta_t)/2 + 1)
psdcall_func = getattr(psd.analytical, psdname)
psd = psdcall_func(length, delta_f, f_min)
detect = detector.Detector(detectorname)

FF_array = func.compute_FF(tb, sg, tau_func, tau_tolerance, psd, nsignal, detect, delta_f, f_min)

fig1, ax1 = plt.subplots()

sg_m1 = [x.m1 for x in sg]
sg_m2 = [x.m2 for x in sg]    
im = plt.scatter(sg_m1, sg_m2, c = FF_array)
#plt.plot(sg_m1, sg_m2, 'x', color='red')
fig1.colorbar(im)
plt.show()
