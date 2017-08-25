#!/usr/bin/env python
import matplotlib.pyplot as plt
from obspy.signal.invsim import corn_freq_2_paz, paz_to_freq_resp
import numpy as np

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

stas = ['FF3','FF6','FF8']


chans =[]
fs =[]
hs=[]
sens=[]
for sta in stas:
    f = open('Results_Swept_Sine' + sta,'r')
    for line in f:
        line= line.split(',')
        fs.append(float(line[4]))
        hs.append(float(line[5]))
        sens.append(float(line[6]))
        chans.append(line[0])
    f.close()


#####################################################Make mean

for triple in zip(fs,hs,sens):
    paz = corn_freq_2_paz(triple[0], triple[1])
    h, f = paz_to_freq_resp(paz['poles'], paz['zeros'], 1., 1./1000., 2**14,True)
    h *= triple[2]
    if 'hmean' not in vars():
        hmean = h
    else:
        hmean +=h

hmean *= 1./9.




fig = plt.figure(1, figsize=(12,12))  
ax = fig.add_subplot(111)
plt.subplots_adjust(hspace=0.001)    
for idx, quad in enumerate(zip(chans, fs, hs, sens)):
    if quad[0] == 'HHZ':
        color='C0'
    elif quad[0] == 'HHN':
        color='C1'
    else:
        color='C2'
    print(quad[1])
    print(quad[2])       
    paz = corn_freq_2_paz(quad[1],quad[2])
    h, f = paz_to_freq_resp(paz['poles'], paz['zeros'], 1., 1./1000., 2**14,True)
    h *= quad[3]
    
    #plt.subplot(2,1,1)
    #if idx < 3:
        #plt.loglog(f, np.abs(h), c=color, label=quad[0])
    #else:
        #plt.loglog(f, np.abs(h), c=color)
    #plt.text(.11, .5*10**5, 'A)', fontsize=24)
    #plt.xlim((.1, 100.))
    #plt.xticks([])
    #plt.legend(loc=4)
    #plt.ylabel('Amplitude (Counts/m/s)')
    #plt.subplot(2,1,2)
    
    #phase = np.unwrap(np.arctan2(-h.imag, h.real))
    #if idx < 3:
        #plt.semilogx(f, phase, c= color, label=quad[0])
    #else:
        #plt.semilogx(f, phase, c= color)
    #plt.text(.11,-.25, 'B)', fontsize=24)
    #plt.legend(loc=4)
    #plt.xlim((.1, 100.))
    #plt.ylabel('Phase (radians)')
    #plt.xlabel('Frequency (Hz)')
    
    plt.subplot(2,1,1)
    if idx <3:
        plt.semilogx(f, 100*(-np.abs(hmean)+np.abs(h))/np.abs(hmean),c=color, label=quad[0])
    else:
        plt.semilogx(f, 100*(-np.abs(hmean)+np.abs(h))/np.abs(hmean),c=color)
    plt.xlim((.1, 100.))
    plt.xticks([])
    plt.text(.11,4.3, 'A)', fontsize=24)
    plt.ylim((-5,5))
    plt.legend(loc=4)
    plt.ylabel('Difference (\%)')
    plt.subplot(2,1,2)
    phasemean = np.unwrap(np.arctan2(-hmean.imag, hmean.real))
    phase = np.unwrap(np.arctan2(-h.imag, h.real))
    if idx < 3:
        plt.semilogx(f, phasemean - phase, c= color)
    else:
        plt.semilogx(f, phasemean - phase, c= color)
    plt.text(.11,.043, 'B)', fontsize=24)
    #plt.legend(loc=4)
    plt.ylim((-.05,.05))
    plt.xlim((.1, 100.))
    plt.ylabel('Phase Difference (radians)')
    plt.xlabel('Frequency (Hz)')
#plt.tight_layout()    
    

    
    
    
    
    
    
plt.savefig('Diffs.jpg',dpi=400)
plt.show()

