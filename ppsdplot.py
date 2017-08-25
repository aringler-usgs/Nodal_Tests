#!/usr/bin/env python
import glob
from obspy.signal import PPSD
import matplotlib.pyplot as plt
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)






files = glob.glob("*.npz")
files = list(files)
files.sort()
fig = plt.figure(1,figsize=(12,12))

for curfile in files:
    ppsd =PPSD.load_npz(curfile)
    #ppsd.plot()
    pers, means = ppsd.get_mean()
    if 'XX' in curfile:
        stalab = 'GS 005 ' + curfile.split('.')[3].replace('PDF','')
    else:
        stalab = curfile.split('.')[0] + ' ' + curfile.split('.')[1] + ' ' + curfile.split('.')[3].replace('PDF','')
    plt.semilogx(1./pers,means, label=stalab,linewidth=2)
 
NHNMper, NHNMpower= get_nhnm()
NLNMper,NLNMpower = get_nlnm()
plt.semilogx(1./NLNMper, NLNMpower, linewidth=3., color='k')
plt.semilogx(1./NHNMper, NHNMpower, linewidth=3., color='k')
plt.xlabel('Frequency (Hz)')
plt.legend(loc=9, ncol=3)
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz$')
plt.gca().invert_xaxis()
plt.xlim((.1,100.))
plt.ylim((-180., -70.))
plt.savefig('MEANPDFnodal.jpg',format='JPEG',dpi=400)
plt.show()
