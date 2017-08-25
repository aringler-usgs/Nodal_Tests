#!/usr/bin/env python
from obspy.core import read, Stream, UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.mlab import csd
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm, PPSD
from obspy.signal.invsim import paz_to_freq_resp
import sys
import glob
from obspy.signal.invsim import corn_freq_2_paz
import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

paz= {'zeros': [0., 0.], 'poles': [-21.99+22.34*1j, -21.99-22.34*1j], 'gain': 0.99813, 'sensitivity': 1.0296*10**9}



#stas = ['YW_4709','YW_4301']
#chans=['DP1','DP2','DPZ']
#for sta in stas:
    #for chan in chans:

        #files = glob.glob('okdata/' + sta +'/*/*/*' + chan + '*.seed')

        #st= Stream()
        #for idx, curfile in enumerate(files):
            #print('On ' + str(idx+1) + ' of ' + str(len(files)))
            #st= read(curfile)
            #if idx == 0:
                #ppsd = PPSD(st[0].stats, paz)
            #ppsd.add(st)
        #ppsd.save_npz(st[0].id + 'PDF.npz')


debug = True
string = 'newdata/output_raw/2017/XX/FF'

pazNEW = corn_freq_2_paz(4.9,.96)
pazNEW['sensitivity'] = 75900.

chans = ['HHZ','HHE', 'HHN']
for idx, chan in enumerate(chans):
    st = Stream()
    for s in range(1,9):

            if debug:
                print('On ' + str(s) + ' ' + chan)
            st += read(string + str(s) + '/' + chan + '.D/XX*')

    stime = UTCDateTime('2017-05-03T01:10:00.0')
    etime = UTCDateTime('2017-05-04T01:10:00.0')
    st.trim(starttime=stime, endtime=etime)
    st.decimate(4)
    st.decimate(2)
    print(st)
    ppsd = PPSD(st[0].stats,pazNEW)
    ppsd.add(st)
    ppsd.save_npz(st[0].id + 'PDF.npz')
