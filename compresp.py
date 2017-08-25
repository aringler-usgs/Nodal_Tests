#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import matplotlib.pyplot as plt
from obspy.signal.invsim import corn_freq_2_paz
from scipy.optimize import fmin_bfgs
import numpy as np
from obspy.signal.cross_correlation import xcorr


import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

paz = {'zeros': [0.j, 0.j, -392. + 0.j, -1960. + 0.j, -1490. + 1740.j, -1490. -1740.j],
    'poles': [-0.03691 + 0.03702j, -0.03691 - 0.03702j, -343. + 0.j, -370. + 467.j, -370. -467.j,
        -836. + 1522.j, -836. -1522.j, -4900. + 4700.j, -4900. - 4700.j, -6900. + 0.j, -15000. + 0.j],
        'gain': 4.344928*10**17, 'sensitivity': 754.3*2.**24/40.}

pazTest= {'zeros': [0., 0.], 'poles': [-21.99+22.34*1j, -21.99-22.34*1j], 'gain': 0.99813, 'sensitivity': 1.0296*10**9}


#sta = 'FF3'
#chans = ['HHZ', 'HHE','HHN']
#stimes = [UTCDateTime('2017-128T19:49:00.0'), UTCDateTime('2017-128T15:19:10'), UTCDateTime('2017-128T15:29:12.')]


#sta = 'FF6'
#chans = ['HHZ','HHE','HHN']
#stimes = [UTCDateTime('2017-128T20:01:00.0'),UTCDateTime('2017-128T15:44:34.0'), UTCDateTime('2017-128T15:53:32.0')]

sta = 'FF8'
chans = ['HHZ','HHE','HHN']
stimes = [UTCDateTime('2017-128T20:14:04.0'),UTCDateTime('2017-128T18:39:20.0'), UTCDateTime('2017-128T19:26:00.0')]





ffile= open('Results_Swept_Sine' + sta,'w')





for idx, pair in enumerate(zip(chans,stimes)):
    stime = pair[1] + 125.
    chan = pair[0]
    #stime = pair[1]
    etime = stime + 225.
    #etime = stime + 300.
    datapath = 'newdata/output_raw/2017/XX/' + sta + '/' + chan + '.D/XX*'
    st = read(datapath)
    st.decimate(5)
    st.decimate(2)
    stRef = read('newdata/Tcmpct/gto.seed')
    
    
    if chan == 'HHE':
       stRef = stRef.select(channel='EH2')
       
    elif chan ==  'HHN':
        stRef = stRef.select(channel='EH1')
        
    else:
        stRef = stRef.select(channel='EH0')
      
    st += stRef
    
    st.sort()
    #st.filter('highpass',freq=1./30.)
    st.trim(stime,etime)

    stRAW = st.copy()
    stRAW.detrend('constant')
    st.taper(0.05)
    if chan == 'HHZ':
        st[0].data *= -1.
    
    # Figure out the sensitivity of the node
    #rat = st[0].std()/st[1].std()
    #print('Here is the ratio of ' + sta + ' ' + chan + ' : ' + str(rat))
    #st[0].data *= 1./rat
    
    
    def resi(x, debug = True):
        
        h=x[1]
        if debug:
            print('Here is h: ' + str(h))
        f=x[0]
        if debug:
            print('Here is f: ' + str(f))
        pazTEMP = corn_freq_2_paz(f,h)
        stTemp = st.copy()
        pazTEMP['sensitivity'] = x[2]
        stTemp[0].simulate(paz_remove=pazTEMP)
        stTemp[1].simulate(paz_remove=paz)
        stTemp.filter('highpass',freq=1./30.)
        comp = sum(np.abs(stTemp[0].data - stTemp[1].data)**2)
        if debug:
            print('Here is the residual: ' + str(comp))
    
        return comp
    
    
    f = 1. /(2.*np.pi / abs(pazTest['poles'][0]))
    h = abs(pazTest['poles'][0].real)/abs(pazTest['poles'][0])
    
    sen = 75900.
    pazTest['sensitivity'] = sen
    
    x = np.array([f, h, sen])
 
    bf = fmin_bfgs(resi,x)
    #else:
    #bf= x
    pazNEW = corn_freq_2_paz(bf[0],bf[1])
    pazNEW['sensitivity'] = bf[2]
    
    stTemp = st.copy()
    
    st[0].simulate(paz_remove=pazNEW)
    st[1].simulate(paz_remove=paz)
    stTemp[0].simulate(paz_remove=pazTest)
    st.filter('highpass',freq=1./30.)
    stTemp.filter('highpass',freq=1./30.)
    rat = st[0].std()/st[1].std()
    print('Here is the ratio of ' + sta + ' ' + chan + ' : ' + str(rat))
    
    ffile.write(chan + ', ' + str(f) + ', ' + str(h) + ', ' + str(sen) + ', ' + str(bf[0]) + ', ' + str(bf[1]) + ', ' + str(bf[2]) + ', ' + str(rat) + ', ' + str(resi(x)) +', ' + str(resi(bf)) + '\n')
    
    index, value = xcorr(st[0],st[1], 300)
    
    print('Here is the index: ' + str(index))
    
    if sta == 'FF3':
        stalab ='005'
    elif sta =='FF6':
        stalab ='006'
    else:
        stalab ='007'
    
    
    fig = plt.figure(1, figsize=(14,14))
    
    plt.subplots_adjust(hspace=0.001)
    t = np.arange(st[0].stats.npts)/st[0].stats.sampling_rate
    plt.subplot(4,1,1)
    plt.title('Swept Sine ' + stalab + ' ' + chan)
    stRAW.normalize()
    plt.plot(t,stRAW[0].data,c='C0')
    plt.xlim((min(t),max(t)))
    plt.ylabel('Counts (Normalized)')
    plt.ylim((-1.2,1.2))
    plt.yticks([-.75,0.,.75])
    plt.xticks([])
    plt.text(5., 1., 'A)')
    plt.subplot(4,1,2)
    plt.plot(t,stRAW[1].data, c='C1')
    plt.ylabel('Counts (Normalized)')
    plt.xlim((min(t),max(t)))
    plt.ylim((-1.2,1.2))
    plt.yticks([-.75, 0., .75])
    plt.text(5., 1., 'B)')
    plt.xticks([])
    plt.subplot(4,1,3)
    
    plt.plot(t,st[1].data*100., label='Trillium Compact', c= 'C1')
    plt.plot(t,stTemp[0].data*100., label= stTemp[0].id + 'Nominal',c='C0')
    plt.plot(t,st[0].data*100.,':',label=st[0].id + ' Perturbed', c='C3')
    plt.text(5., 0.5, 'C)')
    plt.ylabel('Velocity (cm/s)')
    plt.ylim((-.7,.7))
    plt.xlim((min(t),max(t)))
    plt.xticks([])
    plt.subplot(4,1,4)
    #plt.plot(t,st[1].data*100., label='Trillium Compact', c= 'C1')
    plt.plot(t,st[1].data*100.-stTemp[0].data*100., label= stTemp[0].id + 'Nominal',c='C0')
    plt.plot(t,st[1].data*100.-st[0].data*100.,':',label=st[0].id + ' Perturbed', c='C3')
    plt.text(5., 0.27, 'D)')
    
    
    plt.ylabel('Difference (cm/s)')
    plt.ylim((-.35,.35))
    plt.xlim((min(t),max(t)))

    plt.xlabel('Time (s)')
    plt.savefig('SweptSine' + sta + chan + '.jpg', format='JPEG', dpi=400)
    #plt.show()
    plt.clf()
    plt.close()
    #plt.show()
    
ffile.close()
    
    
