import pycbc
from pycbc import waveform, conversions, filter, types, distributions, detector, psd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
import mass as mass
import h5py
import time
import scipy

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
                
def save_matches_HDF(match, filename):
    with h5py.File(filename, 'w') as f:
            f.create_dataset('match', data=match)
    f.close()
    
def read_matches_HDF(filename):
    f = h5py.File(filename, 'r')
    match = np.array(f['match'])
    return match
    
def read_tb(filename, f_min):    
    tb_file = np.loadtxt(filename)
    tb = []
    for i in range(len(tb_file)):
        temp_tau0 = conversions.pycbc.conversions.tau0_from_mass1_mass2(tb_file[i][0], tb_file[i][4], f_min)
        temp_tau3 = conversions.pycbc.conversions.tau3_from_mass1_mass2(tb_file[i][0], tb_file[i][4], f_min)
        temp_mc = conversions.mchirp_from_mass1_mass2(tb_file[i][0], tb_file[i][4])
        temp_q = tb_file[i][0]/tb_file[i][4]
        temp_obj = tb_params(tb_file[i][0], tb_file[i][4], tb_file[i][3], tb_file[i][7], temp_tau0, temp_tau3, temp_mc, temp_q)
        tb.append(temp_obj)
    return tb   

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

def get_mass_ranges_from_TB(tb):
    m1_min = min(x.m1 for x in tb)
    m1_max = max(x.m1 for x in tb)
    m2_min = min(x.m2 for x in tb)
    m2_max = max(x.m2 for x in tb)
    print(m1_min, m1_max, m2_min, m2_max)
    return m1_min, m1_max, m2_min, m2_max
     
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

def generate_template(tb, delta_f, f_min):
    hp, hc = waveform.get_fd_waveform(approximant = "IMRPhenomD",
                                      mass1 = tb.m1,
                                      mass2 = tb.m2,
                                      spin1z = tb.s1z,
                                      spin2z = tb.s2z,
                                      delta_f = delta_f,
                                      f_lower = f_min)
    return hp, hc

def generate_signal(sg, delta_f, f_min): 
    hp, hc = waveform.get_fd_waveform(approximant = "IMRPhenomD",
                                      mass1 = sg.m1,
                                      mass2 = sg.m2,
                                      spin1z = sg.s1z,
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

def resize_wfs(s_f, hp, hc):
    if (len(s_f) > len(hp)):
        hp.resize(len(s_f))
        hc.resize(len(s_f))
    else:
        s_f.resize(len(hp))
    return s_f, hp, hc    


def compute_match(sg, tb, PSD, temp_indices, delta_f, f_min, detect):
    match = []
    start = time.time()
    for i in range(len(temp_indices)):
        ind = temp_indices[i]
        hp_template, hc_template = generate_template(tb[ind], delta_f, f_min)
        hp_signal, hc_signal = generate_signal(sg, delta_f, f_min)  #Loss in match when signal is evaluated outisde the loop
        temp_sf = signal_to_detector_frame(detect, hp_signal, hc_signal, sg)
        
        s_f = temp_sf
        s_f.resize(len(hp_template))
        #s_f, hp, hc = resize_wfs(s_f, hp_template, hc_template)
        match.append(filter.match(hp_template, s_f, psd = PSD, low_frequency_cutoff = f_min)[0])
        #if (match[i] > 0.95):
        #    print(tb[i].tau0, match[i], i)
    end = time.time()
    #print('Exec time', end - start)
    return match   

def fit_tau_envelope(bin_edges, statistic, tau_tolerance):
    #check for NaNs
    idx = np.isfinite(statistic)
    
    f = interp1d(bin_edges[idx], statistic[idx], fill_value="extrapolate")
    
    x_new = np.linspace(min(bin_edges), max(bin_edges), 50)
    y_new = f(x_new) + tau_tolerance
    
    #plt.plot(x_new, y_new, '--', color='orange')
    #plt.show()
    return f


def compute_tauThreshold_envelope(sg_tau0, tau_diff, nbins):
    tau_tolerance = 0.07
    bins = np.linspace(min(sg_tau0), max(sg_tau0), nbins)
    #tau_indices_inbin = np.digitize(sg_tau0, bins)

    statistic, bin_edges, _ = scipy.stats.binned_statistic(sg_tau0, tau_diff, 'max', bins = nbins)
    bin_edges = bin_edges[:-1]
    np.savetxt('injections/tau_threshold.txt', np.c_[bin_edges, statistic])
    tau_func = fit_tau_envelope(bin_edges, statistic, tau_tolerance)
    return tau_func
    
def compute_FF(tb, sg, tau_func, tau_tolerance, psd, nsignal, detect, delta_f, f_min):
    FF_array = []
    recovered_tau = []
    for n in range(nsignal):
        tau0_threshold = tau_func(sg[n].tau0) + tau_tolerance
        template_indices = check_tau0_for_template_generation(tb, sg[n],tau0_threshold)
        if (n%100==0):
            print(n, len(tb), len(template_indices))
        
        match = compute_match(sg[n], tb, psd, template_indices, delta_f, f_min, detect)
        FF_array.append(max(match))
        
    return FF_array