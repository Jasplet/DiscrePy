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
    """
    Function to read sac files based on an input of expected file names (minus the compentend extension).
    It is expected that events will have a North, East and Vertical (BHN,BHE,BHZ) component. If there are not 3 files present then the event is treated as not being read properly.
    Events that are not read properly or that return an exception return a value of False to indicate that there is not sufficient data.
    """
    try:
        st = obspy.core.stream.read(st_id) # Reads all traces that match st_id
        if len(st) == 3: #Check is there are 3 traces in the stream (East,North and Vertical)
            return st
        else:
            return False
    except Exception:
        return False
def save_sac(st,qual,date,time,wbeg,wend):
    """
    Function to trim sac traces once they have been windowed and save these windowed traces for re-use. This is to add easy repeatability
    Arguements:
    st - the stream to process and save
    date - datestamp in the format yyyjjj (julian day)
    time - timestamp in format hhmmss
    wbeg - start of the window
    wend - end of the window
    qual - quaulitativr estimate of seismogram quilty (how clear the SKS arrival is)
    """
    for tr in st:
        ch = tr.stats.channel
        stat = tr.stats.station
        t0 = tr.stats.starttime
        path = '/Users/ja17375/Scripts/Python/Splitting_Codes/SKS_Splitting/Data/Proccessed_Streams'
        tr2 = tr.trim(t0 + wbeg, t0 + wend)
        tr2.write('{}/{}_{}_{:07d}_{:06d}_{}.sac'.format(path,stat,qual,date,time,ch), format='SAC')

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
    rel_SKS = SKS_UTC - st[0].stats.starttime
    return sw.Pair(st[0].data,st[1].data,delta = st[0].stats.delta), rel_SKS

def eigen_plot(eign,fig,**kwargs):
    """
    Function broadly copied from SplitwavePy for plotting splitting measurements from an eigm file.
    The difference being that this verison also returns the figure handle so the figure can be manipulated to add an interactive quality/re-do picker
    """
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
    eigen_plot(split,fig) #Call function eigen_plot to plot splititng messuremnt

    cid = fig.canvas.mpl_connect('key_press_event',interact)
    plt.show(fig)
    fig.canvas.mpl_disconnect(cid)
    return split, pair.wbeg(), pair.wend()

def null():
    """
    Functon to confirm the observation of a null event.
    """
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
    """
    Function for quality assignment.
    Responses differe depending on kep press

    """
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
    else:
        print('Invalid key_press_event, please press a,b,c,n,r or x')

def split_match(date,time,station):
    """
    Function to find and extract the measuremnt from Jack Walpoles splititng data for the same event.
    This matching is done by finding an entry with the same date stamp. Inital testing has shown this to be a unique identifier.
    station MUST be a string of a station code
    date MUST be a int/float of with the format yyyyjjj where j is julian day
    """
#   -------
#   First we need to read in the splitting observations made by Jack Walpole.
#   We also slice out splitting observations just from the station of interest and then reset the indexing so WL_split's indicies start from [0]
#   -------
    raw = pd.read_csv("./Data/Jacks_SKS_RAW.txt",delim_whitespace=True)
    WL_split = raw[(raw['STAT'] == station) & (raw['AUTOQC'] =="split") ]
    WL_split = WL_split.reset_index()
    del WL_split['index']
#   -------
#   Using a Pandas DataFrame we can slice out any rows that match out dare stamp
#   The the iloc function is used to extract the requisite values (Here this is trivial as match should be a single row dataframe, but the values still need to be extracted this way)
#
    match = WL_split[(WL_split['DATE'] == date)] # slicies rows in WL_split that have the same datestamp as date. In theory this should return a single row DataFrame
    if len(match) == 1:
        (fast,dfast,tlag,dtlag,wbeg,wend) = (match.iloc[0]['FAST'],match.iloc[0]['DFAST'],match.iloc[0]['TLAG'],match.iloc[0]['DTLAG'],match.iloc[0]['WBEG'],match.iloc[0]['WEND'])

    elif len(match) == 0:
        print("The provided datestamp does not match any obervations made by JW")
        (fast,dfast,tlag,dtlag,wbeg,wend) = ('NaN','NaN','NaN','NaN','NaN','NaN')
    else:

        print("There has been more than one match, now testing by timestamp also!\n")
        time_test = int(str(time).zfill(6)[0:4]) #Jacks timestamps are only hhmm so I need to strip off the seconds from my timestamps. WARNING it is possible my timestamps are different to Jacks!!
        print('My timestamp {}, Jacks timestamp {}'.format(time_test,match.iloc[0]['TIME']))
        match2 = WL_split[(WL_split['DATE'] == date) & (WL_split['TIME'] == time_test)]
        # print(match2)
        (fast,dfast,tlag,dtlag,wbeg,wend) = (match.iloc[0]['FAST'],match.iloc[0]['DFAST'],match.iloc[0]['TLAG'],match.iloc[0]['DTLAG'],match.iloc[0]['WBEG'],match.iloc[0]['WEND'])

        if len(match2) == 0: #If there is still no match
            (fast,dfast,tlag,dtlag,wbeg,wend) = ('NaN','NaN','NaN','NaN','NaN','NaN')
            print("No match found")
    return fast,dfast,tlag,dtlag,wbeg,wend

###################################################


output_file = open('NEW_Splitting.txt','w+')
output_file.write('STAT DATE TIME STLA STLO EVLA EVLO EVDP GCARC BAZ WBEG WEND FAST DFAST TLAG DTLAG WL_FAST WL_DFAST WL_TLAG WL_DTLAG WL_WBEG WL_WEND QUAL\n')
st_id = []
with open('NEW_read_stream.txt','r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
    for line in reader.readlines():
        line.strip('\n')
        st_id = line[0:-7]
        st = read_sac(st_id+'*')
        # Intialise some global variables which I need to pass things between fucntions (this is probably not be best way to do things but it works!)
        global pair_glob
        global quality
        if st != False: #i.e. if the stream is sufficiently populated and has been read.
            SKS_UTC, t0 = model_SKS(st[0])
            quality = [] # variable to hold Callback key entries for estimated quality of splitting measurements
            date,time = int(str(t0.year)+str(t0.julday).zfill(3)),int(str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+str(t0.second).zfill(2)) #creates time and date stamps
            pair,rel_SKS = st_prep(st = st,trim = 100, f_min = 0.01,f_max = 0.5, SKS = SKS_UTC)
            pair_glob = pair
            # print(' Predicted SKS arrival is at {:5.3f}\n'.format(rel_SKS))
            split, wbeg, wend = measure(pair)
            print('SAC Filename is {}.Window Starts at {:5.2f} and ends at {:5.2f}.\n'.format(line,wbeg,wend,))
#           --------------
#           Now lets find what splitting Jack Walpole measured for this event
#           --------------
            (wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend) = split_match(date,time,"NEW")

            # if quality is not ('x'): #If the quality attribute is not bad (indicated by x)
            filename = '{}_{:07d}_{:06d}.eigm'.format('./Eigm_Files/NEW',date,time)
            attrib = ['stla','stlo','evla','evlo','evdp','gcarc','baz'] #SAC attribute values that I want to extract and print later
            split.save(filename) # Saves splitting measurements
            meas = [wbeg, wend, split.fast, split.dfast, split.lag, split.dlag,wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend ]
            stats = [st[0].stats.sac[i] for i in attrib]
            org = ['NEW',date,time]

            output_file.write('{} {:07d} {:06d} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:06.2f} {:06.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:5.3f} {:4.2f} {} {} {} {} {} {} {}\n'.format(*org,*stats,*meas,quality[0]))
            save_sac(st,quality[0],date,time,wbeg,wend)
            # else:
            #     meas = ['NaN','NaN','NaN','NaN','NaN','NaN']
            #     stats = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN']
            #     org = ['NEW',date, time]

                # output_file.write('{} {:07.0d} {:06.0d} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n'.format(*org,*stats,*meas,quality[0])
        else:
            meas = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN']
            stats = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN']
            org = ['NEW','NaN','NaN']
            quality = ['NaN']
            print('No stream for event',line[0:-7])
            output_file.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22} {23} {24} {25} {26}\n'.format(*org,*stats,*meas,quality[0]))


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
