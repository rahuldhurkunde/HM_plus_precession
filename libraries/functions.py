import pycbc
from pycbc import waveform, conversions, filter, types, distributions, detector, psd, io
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import h5py
import time
import scipy
import multiprocessing
import sys

class tb_params:
        def __init__(self, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, tau0, tau3, mc, q):
                self.m1 = m1
                self.m2 = m2
                self.s1x = s1x
                self.s1y = s1y
                self.s1z = s1z
                self.s2x = s2x
                self.s2y = s2y
                self.s2z = s2z
                self.tau0 = tau0
                self.tau3 = tau3
                self.mc = mc
                self.q = q

class sg_params:
        def __init__(self, m1, m2, s1x, s1y, s1z, s2x, s2y, s2z, tau0, tau3, dist, inc, polarization, right_asc, dec):
                self.m1 = m1
                self.m2 = m2
                self.s1x = s1x
                self.s1y = s1y
                self.s1z = s1z
                self.s2x = s2x
                self.s2y = s2y
                self.s2z = s2z
                self.tau0 = tau0
                self.tau3 = tau3
                self.dist = dist
                self.inc = inc
                self.polarization = polarization
                self.right_asc = right_asc
                self.dec = dec
                
def save_matches_HDF(match, t_ind, filename):
    with h5py.File(filename, 'w') as f:
            f.create_dataset('match', data=match) 
            f.create_dataset('t_ind', data=t_ind)
    f.close()

def read_tb(filename, f_min):    
    tf = h5py.File(filename, 'r')
    print('Template bank --', filename)
    mass1 = tf['mass1']
    mass2 = tf['mass2']
    spin1x = tf['spin1x']
    spin1y = tf['spin1y']
    spin1z = tf['spin1z']
    spin2x = tf['spin2x']
    spin2y = tf['spin2y']
    spin2z = tf['spin2z']
    tb = []
    for i in range(len(mass1)):
        temp_tau0 = conversions.tau0_from_mass1_mass2(mass1[i], mass2[i], f_min)
        temp_tau3 = conversions.tau3_from_mass1_mass2(mass1[i], mass2[i], f_min)
        temp_mc = conversions.mchirp_from_mass1_mass2(mass1[i], mass2[i])
        temp_q = mass1[i]/mass2[i]
        temp_obj = tb_params(mass1[i], mass2[i],
                             spin1x[i], spin1y[i], spin1z[i],
                             spin2x[i], spin2y[i], spin2z[i],
                             temp_tau0, temp_tau3, temp_mc, temp_q)
        tb.append(temp_obj)
    return tb   


def generate_template(tb, delta_f, f_min, approximant_tb, **kwargs):
    hp, hc = waveform.get_fd_waveform(approximant = approximant_tb,
                                      mass1 = tb.m1,
                                      mass2 = tb.m2,
                                      spin1x = tb.s1x,
                                      spin1y = tb.s1y,
                                      spin1z = tb.s1z,
                                      spin2x = tb.s2x,
                                      spin2y = tb.s2y,
                                      spin2z = tb.s2z,
                                      delta_f = delta_f,
                                      f_lower = f_min,
                                      mode_array = kwargs.get('modes'))
    return hp, hc

def generate_signal(sg, delta_f, f_min, approximant_sg, **kwargs):
    hp, hc = waveform.get_fd_waveform(approximant = approximant_sg,
                                      mass1 = sg.m1,
                                      mass2 = sg.m2,
                                      spin1x = sg.s1x,
                                      spin1y = sg.s1y,
                                      spin1z = sg.s1z,
                                      spin2x = sg.s2x,
                                      spin2y = sg.s2y,
                                      spin2z = sg.s2z,
                                      distance = sg.dist,
                                      inclination = sg.inc,
                                      delta_f = delta_f,
                                      f_lower = f_min)
    return hp, hc
            
            
def signal_to_detector_frame(detector, hp_signal, hc_signal, sg):
    fx, fy = detector.antenna_pattern(sg.right_asc, sg.dec, sg.polarization, 1126259462.0)
    signal_detector = hp_signal*fx + hc_signal*fy
    return signal_detector

def check_tau0_for_template_generation(tb, signal, tau0_threshold):
    temp_indices = []
    for i in range(len(tb)):
        if ( abs(tb[i].tau0 - signal.tau0) < tau0_threshold):
            temp_indices.append(i)
    return temp_indices

def real_imag_dot_product(sg, PSD, f_min, delta_f, approximant, nsignal):
	dot_products = []
	for n in range(len(sg)):
		hp, hc = generate_signal(sg[n], delta_f, f_min, approximant)
		hp.resize(len(PSD))
		hc.resize(len(PSD))
		temp = filter.matchedfilter.overlap_cplx(hp, hc, psd = PSD, low_frequency_cutoff = f_min, normalized=True)
		dot_products.append(np.abs(temp.imag))
		if (n % 200 == 0):
			print(n)
	return dot_products    

def flat_tau_envelope(taudiff, taumax, nbins):
    array_len = 200
    taus = np.linspace(0.0, taumax, array_len)
    tau_diff_array = np.full(shape = array_len, fill_value = taudiff)
    filename = 'injections/flat_tau_crawl_%s.txt' %taudiff
    #np.savetxt(filename, np.c_[taus, tau_diff_array])

def compute_tauThreshold_envelope(sg_tau0, tau_diff, nbins):
    bins = np.linspace(min(sg_tau0), max(sg_tau0), nbins)
    statistic, bin_edges, _ = scipy.stats.binned_statistic(sg_tau0, tau_diff, 'max', bins = nbins)
    bin_edges = bin_edges[:-1]
    return bin_edges, statistic
    #np.savetxt('injections/tau_threshold.txt', np.c_[bin_edges, statistic])

def fit_tau_envelope(bin_edges, statistic, tau_tolerance):
    interp_points  = 50     # HARD CODED - try to remove  
    #check for NaNs
    idx = np.isfinite(statistic)
    f = interp1d(bin_edges[idx], statistic[idx], fill_value="extrapolate")
    
    x_new = np.linspace(min(bin_edges), max(bin_edges), interp_points)
    y_new = f(x_new) + tau_tolerance

    #plt.plot(x_new, y_new, '--', color='orange')
    return f

def compute_match(tb, sg, PSD, temp_indices, delta_f, f_min, detect, approximant_tb, approximant_sg, HMs):
    if HMs == True:
        sp, sc = generate_signal(sg, delta_f, f_min, approximant_sg)
    else:
        #print('Only using dominant modes', modes) 
        modes = [[2,2], [2,-2]]
        sp, sc = generate_signal(sg, delta_f, f_min, approximant_sg, modes=modes)  
    s_f = signal_to_detector_frame(detect, sp, sc, sg)
    s_f.resize(len(PSD))
    
    matches = []
    for i in range(len(temp_indices)):  
        ind = temp_indices[i]
        
        hp, hc = generate_template(tb[ind], delta_f, f_min, approximant_tb)
        hp.resize(len(PSD))
        
        temp_match = filter.match(hp, s_f, psd = PSD, low_frequency_cutoff = f_min)[0]
        matches.append(temp_match)
        #if (i%100==0):
            #print(i, temp_match, len(matches))
    return matches   

def compute_FF(tb, sg, tau_func, tau_tolerance, psd, nsignal, detect, delta_f, f_min, approximant_tb, approximant_sg, HMs):
    FF_array = []
    recovered_tau = []
    template_indices = []    

    start = time.time()
    for n in range(nsignal):
        tau0_threshold = tau_func(sg[n].tau0) + tau_tolerance
        template_indices = check_tau0_for_template_generation(tb, sg[n],tau0_threshold)
        if (len(template_indices) == 0):
            sys.exit(' ""Injection outside the parameter region ""')
        else:
            match = compute_match(tb, sg[n], psd, template_indices, delta_f, f_min, detect, approximant_tb, approximant_sg, HMs)
            FF_array.append(max(match))
            recovered_tau.append(tb[template_indices[np.argmax(match)]].tau0)
        #if (n%1==0):
            #end = time.time()
            #print(n, tau0_threshold,len(tb), 'matching templates', len(template_indices), 'Time', end-start)
            #start = time.time()
    return np.array(FF_array), np.array(recovered_tau)

















def resize_wfs(s_f, hp, hc):
    if (len(s_f) > len(hp)):
        hp.resize(len(s_f))
        hc.resize(len(s_f))
    else:
        s_f.resize(len(hp))
    return s_f, hp, hc  


def generate_and_save_templates(tb, delta_f, f_min, approximant_tb, path):
    for n in range(len(tb)):
        hp, hc = generate_template(tb[n], delta_f, f_min, approximant_tb)
        t_ind = np.full(shape = len(hp), fill_value = n)
        filename = '%s/template_%s.hdf' %(path, n)
        try:
            save_template(hp, hc, t_ind, filename)
        except IOError:
            print('path/%s does not exists' %approximant)
            sys.exit()


def save_template(hp, hc, t_ind, filename):
    with h5py.File(filename, 'w') as f:
            f.create_dataset('hp', data=hp) 
            f.create_dataset('hc', data=hc)
            f.create_dataset('t_ind', data=t_ind)
    f.close()
    
def read_template(t_ind, path, delta_f):
    filename = "%s/template_%s.hdf" %(path, t_ind)
    f = h5py.File(filename, 'r')
    hp = types.frequencyseries.FrequencySeries(f['hp'], delta_f=delta_f)
    hc = types.frequencyseries.FrequencySeries(f['hc'], delta_f=delta_f)
    return hp, hc

def mass_samples_from_mc_q(mc_min, mc_max, q_min, q_max, nsignal):
    mc_distribution = mass.MchirpfromUniformMass1Mass2(mc=(mc_min, mc_max-30))   
    q_distribution = mass.QfromUniformMass1Mass2(q=(q_min,q_max-7))
        
    mc_samples = mc_distribution.rvs(size=nsignal)
    q_samples = q_distribution.rvs(size=nsignal)
    m1 = conversions.mass1_from_mchirp_q(mc_samples['mc'],q_samples['q'])
    m2 = conversions.mass2_from_mchirp_q(mc_samples['mc'],q_samples['q'])
    return m1, m2

def mass_samples_from_m1_m2(m1_min, m1_max, m2_min, m2_max, nsignal):
    m1_dist = distributions.Uniform(mass1=(m1_min, m1_max))
    m1 = m1_dist.rvs(size = nsignal)
    
    m2_dist = distributions.Uniform(mass2=(m2_min, m2_max))
    m2 = m2_dist.rvs(size = nsignal)
    return m1, m2

     
def random_params_from_tb(tb, f_min, nsignal):
    mc_min = min(x.mc for x in tb)
    mc_max = max(x.mc for x in tb)
    q_min = min(x.q for x in tb)
    q_max = max(x.q for x in tb)
    m1_min = min(x.m1 for x in tb)
    m1_max = max(x.m1 for x in tb)
    m2_min = min(x.m2 for x in tb)
    m2_max = max(x.m2 for x in tb)
    print( m1_min, m1_max, m2_min, m2_max)
    
    m1, m2 = mass_samples_from_m1_m2(m1_min, m1_max, m2_min, m2_max, nsignal)
    #m1, m2 = mass_samples_from_mc_q(mc_min, mc_max, q_min, q_max, nsignal)
    #tau0 = conversions.tau0_from_mass1_mass2(m1 ,m2, f_min)
    #stau3 = conversions.tau3_from_mass1_mass2(m1, m2, f_min)

    tau0 = []
    tau3 = []
    for k in range(nsignal):
        tau0.append(conversions.tau0_from_mass1_mass2(m1[k][0] ,m2[k][0], f_min))
        tau3.append(conversions.tau3_from_mass1_mass2(m1[k][0] ,m2[k][0], f_min))

    s1z = np.zeros(nsignal)
    s2z = np.zeros(nsignal)
    dist = 200
    inc = 0
    polarization = np.random.uniform(0, 2*np.pi, nsignal)
    right_asc = np.random.uniform(-1, 1, nsignal)
    dec = np.random.uniform(-1, 1, nsignal)
    
    sg = []
    for i in range(nsignal):
        temp_obj = sg_params(m1[i][0], m2[i][0], s1z[i], s2z[i], tau0[i], tau3[i], dist, inc, polarization[i], right_asc[i], dec[i])
        sg.append(temp_obj)
        
    tb_m1 = [x.m1 for x in tb]
    tb_m2 = [x.m2 for x in tb]
    plt.plot(tb_m1, tb_m2, '.', label='templates')
    sg_m1 = [x.m1 for x in sg]
    sg_m2 = [x.m2 for x in sg]
    plt.plot(sg_m1, sg_m2, 'x', color='red', label = 'signals')
    plt.legend()
    plt.show()
    return sg
