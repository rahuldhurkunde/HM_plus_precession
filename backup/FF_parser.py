#!/usr/bin/env python
import sys, os, logging
import argparse as ap
import configparser
import argcomplete
import functions as func
import injections as inj
import pycbc
from pycbc import waveform, conversions, filter, types, distributions, detector, psd
import matplotlib.pyplot as plt
import numpy as np
import scipy
import seaborn as sns
import importlib
import time
import matplotlib.gridspec as gridspec
import Pegasus.DAX3 as dax
import pycbc.workflow as wf
import pycbc.workflow.pegasus_workflow as wdax
from pycbc.workflow import WorkflowConfigParser
import h5py

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

def to_file(path, ifo=None):
    """ Takes a str and returns a pycbc.workflow.pegasus_workflow.File
    instance.
    """
    fil = wdax.File(os.path.basename(path))
    fil.ifo = ifo
    path = os.path.abspath(path)
    fil.PFN(path, "local")
    return fil



# Command line parser
parser = ap.ArgumentParser()
#parser.add_argument('-c', '--config', required=True, help ='specify config file path')
parser.add_argument('-i', '--inference', default = "0", type=int, help ='specify 1 to generate injections (default 0)')
parser.add_argument('-tb', '--template_bank', default = "../banks/vanilla_BBHbank", help ='specify template bank')
parser.add_argument('-inj', '--injections', default = "../injections/new_injections.hdf5", help ='specify injection file')
parser.add_argument('--tau_threshold', default = "../injections/tau_threshold.txt", help ='specify tau thresholds file')
parser.add_argument("--output_dir", default="output/", help="Path to output directory.")

# add option groups
wf.add_workflow_command_line_group(parser)

# Parser command line
args, remaining_args = parser.parse_known_args()
if args.config_files[0]:
	config = configparser.ConfigParser()
	config.read([args.config_files[0]])

f_min, delta_f, delta_t, psdname, detectorname, nsignal, nbins_tau, tau_tolerance = determine_values_from_parsing_config(config)
pycbc.init_logging(True)

# create workflow and sub-workflows
workflow = wf.Workflow(args, "gw")

# make data output and results directories
wf.makedir(args.output_dir)

infhand = to_file(args.config_files[0])

if args.inference == 1:
	#construct executable for pycbc_generate_injections
	injections_exe = wf.Executable(workflow.cp, "injections", ifos=workflow.ifos,
                              out_dir=args.output_dir)

	path = args.config_files[1]
	inj_config = to_file(path)

	#make node for running injection-generator
	node = injections_exe.create_node()
	node.add_input_opt('--config-file', inj_config)
	inj_file = node.new_output_file_opt(workflow.analysis_time, ".hdf", "--output-file")
	workflow += node

# write dax
workflow.save('gw.dax')
exit()

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
