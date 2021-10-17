#Script to read injections using an HDF5 file. Also, contains functions to generate random parameters within a given template bank region

import pycbc
from pycbc import waveform, conversions, filter, types, distributions, detector, psd
import numpy as np
import mass as mass
import functions as func
import h5py

class tb_params:
        def __init__(self, m1, m2, s1z, s2z, tau0, tau3, mc, q):
                self.m1 = m1
                self.m2 = m2
                self.s1z = s1z
                self.s2z = s2z
                self.tau0 = tau0
                self.tau3 = tau3
                self.mc = mc
                self.q = q

class sg_params:
        def __init__(self, m1, m2, s1z, s2z, tau0, tau3, dist, inc, polarization, right_asc, dec):
                self.m1 = m1
                self.m2 = m2
                self.s1z = s1z
                self.s2z = s2z
                self.tau0 = tau0
                self.tau3 = tau3
                self.dist = dist
                self.inc = inc
                self.polarization = polarization
                self.right_asc = right_asc
                self.dec = dec

def write_injections_HDF(sg, filename):
    m1 = np.array([x.m1 for x in sg])
    m2 = np.array([x.m2 for x in sg])
    s1z = np.array([x.s1z for x in sg])
    s2z = np.array([x.s2z for x in sg])
    tau0 = np.array([x.tau0 for x in sg])
    tau3 = np.array([x.tau3 for x in sg])
    dist = np.array([x.dist for x in sg])
    inc = np.array([x.inc for x in sg])
    polarization = np.array([x.polarization for x in sg])
    right_asc = np.array([x.right_asc for x in sg])
    dec = np.array([x.dec for x in sg])

    with h5py.File(filename, 'w') as f:
            f.create_dataset("m1", data=m1)
            f.create_dataset("m2", data=m2)
            f.create_dataset("s1z", data=s1z)
            f.create_dataset("s2z", data=s2z)
            f.create_dataset("tau0", data=tau0)
            f.create_dataset("tau3", data=tau3)
            f.create_dataset("dist", data=dist)
            f.create_dataset("inc", data=inc)
            f.create_dataset("polarization", data=polarization)
            f.create_dataset("right_asc", data=right_asc)
            f.create_dataset("dec", data=dec)
    f.close()
    
def read_injections_HDF(filename, nsignal):
    hf = h5py.File(filename, 'r')
    sg = []
    for i in range(nsignal):
        temp_obj = sg_params(hf['m1'][i], 
                            hf['m2'][i], 
                            hf['s1z'][i], 
                            hf['s2z'][i], 
                            hf['tau0'][i], 
                            hf['tau3'][i], 
                            hf['dist'][i], 
                            hf['inc'][i], 
                            hf['polarization'][i], 
                            hf['right_asc'][i], 
                            hf['dec'][i])
        sg.append(temp_obj)
    return sg


def generate_injections(delta_f, delta_t, f_min):
    tb = func.read_tb("banks/vanilla_BBHbank", f_min)
    nsignal = 1000
    sg = func.random_params_from_tb(tb, f_min, nsignal)
    filename = 'random_injections.hdf5'
    #write_injections_HDF(sg, filename)


