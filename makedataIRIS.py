#!/usr/bin/env python
from obspy.core import read, Stream
import glob
import sys

files = glob.glob('newdata/output_raw/2017/XX/*/*/XX*')


for curfile in files:
    print(curfile)
    st=read(curfile)
    for tr in st:
        tr.stats.network = 'GS'
        if tr.stats.station in ['FF1', 'FF2', 'FF3']:
            tr.stats.station ='005'
        elif tr.stats.station in ['FF4','FF5','FF6']:
            tr.stats.station ='006'
        elif tr.stats.station in ['FF7','FF8']:
            tr.stats.station = '007'
        tr.stats.channel = 'E' + tr.stats.channel[1:]
    st.write(curfile + '.mseed' ,format='MSEED')
            

