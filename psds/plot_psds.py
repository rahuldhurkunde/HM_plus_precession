import numpy as np
import matplotlib.pyplot as plt
import sys
import pycbc
from pycbc import psd


#nplots = len(sys.argv) - 1
#print("No. of psds given", nplots)
#
#if (nplots%2 == 0):
#	fig, axs = plt.subplots(int(nplots/2), 2)
#else:
#	fig, axs = plt.subplots(int((nplots+1)/2), 2)
#
#for n in range(1, nplots+1):
#	filename = sys.argv[n]
#	freq = np.loadtxt(filename)[:,0]
#	asd = np.loadtxt(filename)[:,1]
#	ax = axs.ravel()[n-1]
#	ax.set_xlim([0, 1000])
#	ax.plot(freq, asd, label=filename)
#	ax.set_title(filename)
#
labels = ['Advanced LIGO design', 'A+', 'LIGO Voyager', 'Cosmic explorer']

fig, ax = plt.subplots()
nplots = len(sys.argv) - 1
for n in range(0, nplots):
	filename = sys.argv[n+1]
	freq = np.loadtxt(filename)[:,0]
	asd = np.loadtxt(filename)[:,1]
	ax.plot(freq, asd, label=labels[n])
	#ax.set_xlim([7, 512])
	#ax.set_ylim([0, 3.2*10**(-20)])


plt.grid()
plt.yscale('log')
plt.xscale('log')
plt.yticks(fontsize=13)
plt.xticks(fontsize=13)

plt.xlabel('Frequency (Hz)', fontsize=15)
plt.ylabel('$\mathdefault{Strain} / \sqrt{\mathdefault{Hz}}$', fontsize=15)
plt.legend()
plt.savefig('PSDs.png', dpi=600, bbox_inches='tight')
plt.show()
