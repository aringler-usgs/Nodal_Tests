#!/usr/bin/env python
from obspy.core import read, Stream, UTCDateTime
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.mlab import csd
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
from obspy.signal.invsim import paz_to_freq_resp
import sys

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

paz= {'zeros': [0., 0.], 'poles': [-21.99+22.34*1j, -21.99-22.34*1j], 'gain': 0.99813, 'sensitivity': 1.0296*10**9}

pazSTS2={'zeros': [0., 0., -15.15, -176.6, -463 -430.5*1j, -463+430.5*1j], 'poles': [-0.037 +0.037*1j, -0.037-0.037*1j, -15.64, -97.34+400.7*1j, 
        -97.34-400.7*1j, -374.8, -520.3, -1053.-1005.*1j,  -1053.+ 1005.*1j, -1330., -255.097], 'gain':3.468*10**17, 'sensitivity': 20000.*(2.**26)/40.}


def computeresp(resp,delta,lenfft):
    respval,freq = paz_to_freq_resp(resp['poles'],resp['zeros'],resp['sensitivity']*resp['gain'],t_samp = delta, 
        nfft=lenfft,freq = True)
    idx = np.argmin(np.abs(freq-.1))

    respval = np.absolute(respval*np.conjugate(respval))
    respval *= 1./respval[idx]
    respval = respval[1:]*resp['sensitivity']**2
    respval = respval.real
    return respval

def cp(data1,data2,lenfft,lenol,delta):
    sr = 1/delta
    cpval,fre = csd(data1,data2,NFFT=lenfft,Fs=sr,noverlap=lenol,scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:]
    return cpval, fre 
    
    
lenfft=4096*2*2*2
lenol=512*2*2*2

def selfnoise(st):
    delta = st[0].stats.delta
    instresp = computeresp(paz, delta, lenfft)
    (p11, f) = cp(st[0],st[0],lenfft,lenol,delta)
    (p22, f) = cp(st[1],st[1],lenfft,lenol,delta)
    (p33, f) = cp(st[2],st[2],lenfft,lenol,delta)

    (p21, f) = cp(st[1],st[0],lenfft,lenol,delta)
    (p13, f) = cp(st[0],st[2],lenfft,lenol,delta)
    (p23, f) = cp(st[1],st[2],lenfft,lenol,delta)

    n11 = (p11 - p21*p13/p23)
    n22 = (p22 - np.conjugate(p23)*p21/np.conjugate(p13))
    n33 = (p33 - p23*np.conjugate(p13)/p21)
    n11 = 10.*np.log10(((2*np.pi*f)**2)*n11/instresp)
    n22 = 10.*np.log10(((2*np.pi*f)**2)*n22/instresp)
    n33 = 10.*np.log10(((2*np.pi*f)**2)*n33/instresp)
    p11= 10.*np.log10(((2*np.pi*f)**2)*p11/instresp)
    p22= 10.*np.log10(((2*np.pi*f)**2)*p22/instresp)
    p33= 10.*np.log10(((2*np.pi*f)**2)*p33/instresp)
    n = [n11, n22, n33]
    p=[p11, p22, p33]
    
    return n, p, f

debug = True
string = 'newdata/output_raw/2017/XX/FF'

# decide if it is high-frequency or not
HF = True

chans = ['HHZ','HHE', 'HHN']
fig = plt.figure(1, figsize=(12,12))
plt.subplots_adjust(hspace=0.001)
for idx, chan in enumerate(chans):
    st = Stream()
    for s in range(1,9):

            if debug:
                print('On ' + str(s) + ' ' + chan)
            st += read(string + str(s) + '/' + chan + '.D/XX*')

    stime = UTCDateTime('2017-05-03T01:10:00.0')
    etime = UTCDateTime('2017-05-03T02:10:00.0')
    st.trim(starttime=stime, endtime=etime)



    nsHF, psHF,fHF = selfnoise(st)
    st.decimate(5)

    st.decimate(2)

    st.decimate(5)

    print(st)
    ns, ps,f = selfnoise(st)




    #st.plot()
    if chan == 'HHE':
        sts2chan = 'BH2'
    elif chan == 'HHN':
        sts2chan ='BH1'
    elif chan =='HHZ':
        sts2chan = 'BH0'
    st2 = read('/msd/XX_TST1/2017/' + str(stime.julday).zfill(3) + '/00_' + sts2chan + '*')
    st2.detrend('constant')

    st2.trim(starttime=stime,endtime=etime)


    pref,f2 = cp(st2[0].data,st2[0].data,lenfft,lenol,st2[0].stats.delta)
    instresp = computeresp(pazSTS2, st2[0].stats.delta, lenfft)

    pref = 10.*np.log10(((2*np.pi*f2)**2)*pref/instresp)






    NLNMper,NLNMpower = get_nlnm()
    NHNMper,NHNMpower = get_nhnm()

    pm = ps[0]
    pm = np.vstack((ps[1], pm))
    pm = np.vstack((ps[2], pm))


    pm= np.mean(pm, axis=0)

    nm = ns[0]
    nm = np.vstack((ns[1], nm))
    nm = np.vstack((ns[2], nm))


    nm= np.mean(nm, axis=0)

    pmHF = psHF[0]
    pmHF = np.vstack((psHF[1], pmHF))
    pmHF = np.vstack((psHF[2], pmHF))


    pmHF= np.mean(pmHF, axis=0)

    nmHF = nsHF[0]
    nmHF = np.vstack((nsHF[1], nmHF))
    nmHF = np.vstack((nsHF[2], nmHF))

    nmHF= np.mean(nmHF, axis=0)

    N= 20
    pmHF=np.convolve(pmHF, np.ones((N,))/N, mode='same')
    fHF=np.convolve(fHF,np.ones((N,))/N, mode='same')
    nmHF=np.convolve(nmHF,np.ones((N,))/N, mode='same')
    pm=np.convolve(pm, np.ones((N,))/N, mode='same')
    f=np.convolve(f,np.ones((N,))/N, mode='same')
    nm=np.convolve(nm,np.ones((N,))/N, mode='same')
    f2=np.convolve(f,np.ones((N,))/N, mode='same')
    pref=np.convolve(pref,np.ones((N,))/N, mode='same')

    minfre = .16
    maxfre = .3


    scale = np.mean(pref[(f2 >= minfre) & (f2 <= maxfre)])- np.mean(ps[0][(f >= minfre) & (f <= maxfre)])
    print(scale)


    scaleHF = np.mean(pm[(f >= 10.) & (f <= 11.)]) - np.mean(pmHF[(fHF >= 10.) & (fHF <= 11.)])




    numb = 20
    
    plt.subplot(3,1,idx+1)
    if HF:
        for idx in range(3):
            plt.semilogx(fHF[numb:-numb],psHF[idx][numb:-numb]+scale+scaleHF,color='g')
            plt.semilogx(fHF[numb:-numb],nsHF[idx][numb:-numb]+scale+scaleHF,color='c')
        plt.semilogx(fHF[numb:-numb],pmHF[numb:-numb]+scale+scaleHF,color='b')
        plt.semilogx(fHF[numb:-numb],nmHF[numb:-numb]+scale+scaleHF,color='r')
        
    else:
        for idx in range(3):
            plt.semilogx(f[numb:-numb],ps[idx][numb:-numb]+scale+scaleHF,color='g')
            plt.semilogx(f[numb:-numb],ns[idx][numb:-numb]+scale+scaleHF,color='c')
        plt.semilogx(f[numb:-numb],pm[numb:-numb]+scale,color='b',label='Mean Nodal Noise')
        plt.semilogx(f[numb:-numb],nm[numb:-numb]+scale,color='r', label='Mean Nodal Incoherent Noise')
        plt.semilogx(f2[numb:-numb],pref[numb:-numb],c='.5',label='STS-2 Reference')
    

    plt.semilogx(1./NLNMper, NLNMpower, linewidth=2., color='k')
    #plt.semilogx(1./NHNMper, NHNMpower, linewidth=2., color='k')
    if idx == 0:
        let = 'A)'
    elif idx == 1:
        let = 'B)'
    else:
        let = 'C)'
    if HF:
        plt.text(5.5, -110, let + '  ' + chan, fontsize=28)
    else:
        plt.text(.11, -110, let + '  ' + chan, fontsize=28)
    if HF:
        plt.xlim((1000.,5.))
    else:
        plt.xlim((30.,.1))
    if idx < 2:
        plt.xticks([])
    
    plt.ylim((-200., -100.))
    plt.yticks([-180., -160., -140., -120.])
    if idx == 1:
        plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz)$')
    plt.xlabel('Frequency (Hz)')
    plt.gca().invert_xaxis()
    #plt.title(chan + ' Self-Noise Estimates for ' + str(stime.year) + ' ' + str(stime.julday).zfill(3) + str(stime.hour).zfill(2) + ':' + str(stime.minute).zfill(2) + ' Duration 1 Hour')
    #plt.legend()
    if HF:
        plt.savefig('Self_noiseHFESUPP.jpg', format='JPEG', dpi=400)
    else:
        plt.savefig('Self_noiseLPESUPP.jpg', format='JPEG', dpi=400)

plt.show()
