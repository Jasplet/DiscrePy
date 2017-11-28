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
model = obspy.taup.tau.TauPyModel(model="iasp91")

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
        work = True
    except Exception:
        print("Exception Encountered for event", i)
        row = str(i)+' N/A N/A N/A N/A N/A N/A N/A N/A N/A N/A N/A \n'
        work = False
#def Split_Measure(no_events):
for i in range(0,101):
    st_id = "NEW_" + str(i).zfill(2) + "_" + "*.sac" # Generate expected file names. Wildcard used to catch all 3 channels
    filename = "./Splitting/NEW_" + str(i).zfill(2) + ".eigm"
    read_sac(st_id)
    if (work == True and len(st) == 3 ): #Traces for event have been succesfully read so lets try to measure splittiing!
        SKS = model.get_travel_times(st[0].stats.sac.evdp,st[0].stats.sac.gcarc,["SKS"])[0].time # Calculates expected SKS arrival for the current event
        t0 = st[0].stats.starttime # Start time of stream for event. This should be the event origin time.        st.filter("bandpass",freqmin=0.01, freqmax= 0.5,corners=2,zerophase=True)
        SKS_arr = obspy.core.utcdatetime.UTCDateTime(st[0].stats.starttime + SKS)# SKS arrival time relative to trace start as a UTCDateTime object
        st.filter("bandpass",freqmin=0.01, freqmax= 0.5,corners=2,zerophase=True) # Zerophase bandpass filter of the streams
        st.trim(starttime = SKS_arr - 120, endtime = SKS_arr +180)
        east = st[0].data
        north = st[1].data
        split_pair = sw.Pair(east,north,delta = st[0].stats.delta) # Creates the Pair object (East compenent,North compenent,sample_interval)
        #split_pair.set_window(SKS_arr - 30,SKS_arr + 30)
        split_pair.plot(pick=True ,marker =  120) # Plots traces and particle motion, with a marker for the SKS arrival. Plot window also allows manual picking of the window.
        split = sw.EigenM(split_pair,lags=(4,) )
        split.save(filename) # Saves splitting measurements
        window1 = split_pair.wbeg()
        window2 = split_pair.wend()
        date = str(t0.year)+" "+ str(t0.month).zfill(2)+" "+str(t0.day).zfill(2) +" " # Format date information so that is inteligible when output
        t = str(t0.hour).zfill(2)+" "+str(t0.minute).zfill(2)+" "+str(t0.second).zfill(2) # Formating time infomation for printing
        row = str(i)+' '+date+' '+t+' NEW '+str(split.fast)+' '+str(split.dfast)+' '+str(split.lag)+' '+str(split.dlag)+ ' '+ str(window1) + ' '+ str(window2) +'\n' #Row of data to be written to output textfile
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
