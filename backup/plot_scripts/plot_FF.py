import pycbc 
import numpy as np
import matplotlib.pyplot as plt
import h5py

f = h5py.File('TF2-FF-0-10.hdf', 'r')
FF = np.array(f['FF'])
m1 = np.array(f['mass1'])
m2 = np.array(f['mass2'])

fig1, ax1 = plt.subplots()

im = plt.scatter(m1, m2, c = FF)
fig1.colorbar(im)
plt.show()

