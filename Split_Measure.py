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

#### Reading downloaded traces and Windowing data
# Generate Travel time model for SKS phase for all events prior to Windowing


# arrival = model.get_travel_times(depth,dist,["SKS"])
# no_N = int(os.popen("ls *BHN.sac | wc -l").read())
# no_E = int(os.popen("ls *BHE.sac | wc -l").read())
# no_Z = int(os.popen("ls *BHZ.sac | wc -l").read())
# num_s = max(no_N,no_E) # Find which channel has the max number of traces downloaded (i.e what is the maximum number of events we can measure splitting for)

output_file = open('NEW_Splitting.txt','w')
output_file.write('ID YEAR MON DAY HOUR MIN SEC STAT FAST DFAST TLAG DTLAG WBEG WEND \n')

def read_sac(st_id):
    try:
        st = obspy.core.stream.read('./Data/' + st_id) # Reads all traces that match st_id
        if len(st) == 3: #Check is there are 3 traces in the stream (East,North and Vertical)
            return True
        else:
            return False
    except Exception:
        return False

def model_SKS(tr):
    model = obspy.taup.tau.TauPyModel(model="iasp91")
    SKS = model.get_travel_times(tr.stats.sac.evdp,tr.stats.sac.gcarc,["SKS"])[0].time
    t0 = st[0].stats.starttime # Start time of stream for event. This should be the event origin time.
    UTC = obspy.core.utcdatetime.UTCDateTime(t0 + SKS)# SKS arrival time relative to trace start as a UTCDateTime object
    return SKS,t0,UTC
    #def Split_Measure(no_events):

def st_prep(st,trim_beg,trim_end,f_min,f_max):
    """
    Prepares Stream for spltting analysis and creates the Pair object
    """
    st.filter("bandpass",freqmin= f_min, freqmax= f_max,corners=2,zerophase=True) # Zerophase bandpass filter of the streams
    return st.trim(starttime = SKS_UTC - trim_beg, endtime = SKS_UTC + trim_end)

def window_trace(pair,ext,*args):
    """
    window_trace(pair,ext,st,f_min,_f_max,SKS)
    Constructs Pair Object for a Given stream and the generates the splitwavepy window Picker
    If Pair does not exists then pass 'pair = None' and st,f_min,f_max
    st = stream , f_min/max are the min/max filter frequencies.
    ext refers to the trim size and also the position of the SKS marker
    Up to 3 arguements accepted. If pair = none then a stream MUST also be provided.
    The SKS marker is hard coded in as 150 seconds from the start of the trimmed trace. This is the default value.

    """
    if pair is None:
        st.filter("bandpass",freqmin= f_min, freqmax= f_max,corners=2,zerophase=True) # Zerophase bandpass filter of the streams
        st.trim(starttime = SKS - trim, endtime = SKS + trim) #Trims Stream cenetered around the predicited SKS arrival
        pair = sw.Pair(st[0].data,st[1].data,delta = st[0].stats.delta) # Creates the Pair object (East compenent,North compenent,sample_interval)

    pair.plot(pick=True,marker = 150) # Plots the window picker.
    return pair

    if len(args) > 6:
        raise Exception('Too Many Arguements Provided')

for i in range(0,101):
    st_id = "NEW_" + str(i).zfill(2) + "_" + "*.sac" # Generate expected file names. Wildcard used to catch all 3 channels
    filename = "./Splitting/NEW_" + str(i).zfill(2) + ".eigm"
    isread = read_sac(st_id)

    if isread == True:
        #Traces for event have been succesfully read so lets try to measure splittiing!
        (SKS, origin, SKS_UTC) = model_SKS(st[0])
        pair = window_trace(pair = None,150,st = st, f_min = 0.01,f_max = 0.5, SKS = SKS_UTC)

        split = sw.EigenM(pair,lags=(4,) )
        split.plot()
        ## Callback key entries for estimated quality of splitting measurements
            #if key is (A,B,C):
            #   split.save(filename) # Saves splitting measurements
                window1 = split_pair.wbeg()
                window2 = split_pair.wend()
                date = str(t0.year)+" "+ str(t0.month).zfill(2)+" "+str(t0.day).zfill(2) +" " # Format date information so that is inteligible when output
                t = str(t0.hour).zfill(2)+" "+str(t0.minute).zfill(2)+" "+str(t0.second).zfill(2) # Formating time infomation for printing
                row = str(i)+' '+date+' '+t+' NEW '+str(split.fast)+' '+str(split.dfast)+' '+str(split.lag)+' '+str(split.dlag)+ ' '+ str(window1) + ' '+ str(window2) +'\n' #Row of data to be written to output textfile
            #elif key is X:
                #window_trace(pair, 150)
                #measure splitting again 

    row = str(i)+' N/A N/A N/A N/A N/A N/A N/A N/A N/A N/A N/A \n'
    output_file.write(row)

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
