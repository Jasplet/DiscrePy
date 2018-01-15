#! /usr/bin/env python
## This is a script which takes downloaded traces and Calculates SKS splitting
import obspy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import obspy.signal.filter as filt
import obspy.core.utcdatetime
import glob
import os.path
import splitwavepy as sw
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
#### Reading downloaded traces and Windowing data
# Generate Travel time model for SKS phase for all events prior to Windowing


# arrival = model.get_travel_times(depth,dist,["SKS"])
# no_N = int(os.popen("ls *BHN.sac | wc -l").read())
# no_E = int(os.popen("ls *BHE.sac | wc -l").read())
# no_Z = int(os.popen("ls *BHZ.sac | wc -l").read())
# num_s = max(no_N,no_E) # Find which channel has the max number of traces downloaded (i.e what is the maximum number of events we can measure splitting for)

def read_sac(st_id):
    try:
        st = obspy.core.stream.read(st_id) # Reads all traces that match st_id
        if len(st) == 3: #Check is there are 3 traces in the stream (East,North and Vertical)
            return st
        else:
            return False
    except Exception:
        return False

def model_SKS(tr):
    """
    Function to run TauP traveltime models for the SKS phase.
    Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object
    """
    model = obspy.taup.tau.TauPyModel(model="iasp91")
    SKS = model.get_travel_times(tr.stats.sac.evdp,tr.stats.sac.gcarc,["SKS"])[0].time
    t0 = st[0].stats.starttime # Start time of stream for event. This should be the event origin time.
    UTC = obspy.core.utcdatetime.UTCDateTime(t0 + SKS)# SKS arrival time relative to trace start as a UTCDateTime object
    return UTC, t0
    #def Split_Measure(no_events):

def st_prep(st,trim,f_min,f_max,SKS):
    """
    Prepares Stream for spltting analysis (by bandpass filtering and trimming) and then creates the Pair object
    """
    st.filter("bandpass",freqmin= f_min, freqmax= f_max,corners=2,zerophase=True) # Zerophase bandpass filter of the streams
    st.trim(starttime = SKS_UTC - trim, endtime = SKS_UTC + trim)
    return sw.Pair(st[0].data,st[1].data,delta = st[0].stats.delta)

def eigen_plot(eign,fig,**kwargs):

    # setup figure and subplots
    #fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(2, 3,width_ratios=[1,1,2])
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[0,1])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[1,1])
    ax4 = plt.subplot(gs[:,2])

    # data to plot
    d1 = eign.data.chop()
    d1f = eign.srcpoldata().chop()
    d2 = eign.data_corr().chop()
    d2s = eign.srcpoldata_corr().chop()

    # flip polarity of slow wave in panel one if opposite to fast
    # d1f.y = d1f.y * np.sign(np.tan(eign.srcpol()-eign.fast))

    # get axis scaling
    lim = np.abs(d2s.data()).max() * 1.1
    ylim = [-lim,lim]

    # original
    d1f._ptr(ax0,ylim=ylim,**kwargs)
    d1._ppm(ax1,lims=ylim,**kwargs)
    # corrected
    d2s._ptr(ax2,ylim=ylim,**kwargs)
    d2._ppm(ax3,lims=ylim,**kwargs)

    # error surface
    if 'vals' not in kwargs:
        # kwargs['vals'] = (eign.lam1 - eign.lam2) / eign.lam2
        # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
        kwargs['vals'] = eign.lam1 / eign.lam2
        kwargs['title'] = r'$\lambda_1 / \lambda_2$'

    # add marker and info box by default
    if 'marker' not in kwargs: kwargs['marker'] = True
    if 'info' not in kwargs: kwargs['info'] = True
    if 'conf95' not in kwargs: kwargs['conf95'] = True
    eign._psurf(ax4,**kwargs)
    # neaten
    plt.tight_layout()

def measure(pair):
    """
    Function for Picking the window for a provded pair object and then measure the splitting
    """
    pair.plot(pick=True,marker = 100) # Plots the window picker.
    split = sw.EigenM(pair,lags=(4,) )

    fig = plt.figure(figsize=(12,6))
    eigen_plot(split,fig)

    cid = fig.canvas.mpl_connect('key_press_event',interact)
    plt.show(fig)
    fig.canvas.mpl_disconnect(cid)
    return split, pair.wbeg(), pair.wend()

def null():
    var = input("Null Measurement (y/n)")
    if var is 'y':
        quality.append('null')
        print("Null Confirmed")
        plt.close()
    elif var is 'n':
        print("Ok, repeating measurement")
        measure(pair_glob)
    else:
        print("Invalid response, try again")
        null()

def interact(event):
    if event.key in ['a','b','c']:
        quality.append(event.key)
        plt.close()
    elif event.key is 'n':
        null()
    elif event.key is 'x':
        quality.append('x')
        plt.close()
    elif event.key == 'r':
        plt.close()
        measure(pair_glob)


###################################################


output_file = open('NEW_Splitting.txt','w')
output_file.write('STAT DATE TIME STLA STLO EVLA EVLO EVDP GCARC BAZ FAST DFAST TLAG DTLAG WBEG WEND QUAL\n')
st_id = []
with open('NEW_read_stream.txt','r') as reader:
    for line in reader.readlines():
        line.strip('\n')
        st_id = line[0:-7]
        st = read_sac(st_id+'*')
        # Intialise some global variables which I need to pass things between fucntions (this is probably not be best way to do things but it works!)
        global pair_glob
        global quality
        if st != False: #i.e. if the stream is sufficiently populated and has been read.
            SKS_UTC, t0 = model_SKS(st[0])
            quality = []
            pair = st_prep(st = st,trim = 100, f_min = 0.01,f_max = 0.5, SKS = SKS_UTC)
            pair_glob = pair
            split, wbeg, wend = measure(pair)
            print(wbeg,wend)
            ## Callback key entries for estimated quality of splitting measurements
            if quality is not ('x'): #If the quality attribute is not bad (indicated by x)
                filename = '{}_{:04d}_{:02d}_{:02d}_{:02d}_{:02d}_{:02d}.eigm'.format('./Splitting/NEW',t0.year,t0.month,t0.day,t0.hour,t0.minute,t0.second)
                attrib = ['stla','stlo','evla','evlo','evdp','gcarc','baz']
                split.stla = st[0].stats.sac[attrib[0]]
                split.stlo = st[0].stats.sac[attrib[1]]
                split.evla = st[0].stats.sac[attrib[2]]
                split.evlo = st[0].stats.sac[attrib[3]]
                split.evdp = st[0].stats.sac[attrib[4]]
                split.gcarc = st[0].stats.sac[attrib[5]]
                split.baz = st[0].stats.sac[attrib[6]]
                split.save(filename) # Saves splitting measurements
                meas = [split.fast, split.dfast, split.lag, split.dlag, wbeg, wend]
                stats = [st[0].stats.sac[i] for i in attrib]
                org = ['NEW',int(str(t0.year)+str(t0.julday).zfill(3)),int(str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+str(t0.second).zfill(2))]
                row = tuple(org + stats + meas) #Row of data to be written to output textfile
                row = row + (quality,)
                print(row)
                output_file.write('{} {:07.0d} {:06.0d} {:06.2f} {:06.2f} {:06.2f} {:06.2f} {:06.3f} {:06.2f} {:06.2f} {:4.2f} {:4.2f} {:5.3f} {:4.2f} {:4.2f} {:4.2f} {}\n'.format(*org,*stats,*meas,quality[0]))
            else:
                meas = ['N/A','N/A','N/A','N/A','N/A','N/A']
                stats = ['N/A','N/A','N/A','N/A','N/A','N/A','N/A']
                org = ['NEW',int(str(t0.year)+str(t0.julday).zfill(3)),int(str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+str(t0.second).zfill(2))]
                row = tuple(org + stats + meas)

                output_file.write('{} {:07.0d} {:06.0d} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n'.format(*out))
        else:
            meas = ['N/A','N/A','N/A','N/A','N/A','N/A']
            stats = ['N/A','N/A','N/A','N/A','N/A','N/A','N/A']
            org = ['NEW','N/A','N/A','N/A','N/A','N/A','N/A']
            row = tuple(org + stats + meas + ['NoStream'])
            print('No stream for event',line[0:-7])
            output_file.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20}\n'.format(*row))


output_file.close()

#data = [split1.FAST.data,split1.TLAG.data,fast.data,tlag.data]
#print(data)
# index = range(0,len(split1.FAST))
# columns = ["Jacks_Fast","Jacks_Lag","My_Fast","My_Lag"]
# df = pd.DataFrame(index=index,columns=columns)
# df = df.fillna(0)
# df.Jacks_Fast = data[0]
# df.Jacks_Lag = data[1]
# df.My_Fast = data[2]
# df.My_Lag = data[3]
# print(df)
# fast_lag = df[(df['My_Fast'] != 0) & (df['My_Lag'] != 0)]
# print(fast_lag)
