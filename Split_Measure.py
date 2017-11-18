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


for i in range(99,100):
    st_id = "NEW_" + str(i).zfill(2) + "_" + "*.sac" # Generate expected file names. Wildcard used to catch all 3 channels

    filename = "NEW_Splitting_" + str(95).zfill(2) + ".eigm"

    try:
        st = obspy.core.stream.read('/Users/ja17375/Scripts/Python/SKS_Splitting/Data/' + st_id) # Reads all traces that match st_id
        work = True

    except Exception:
        print("Exception Encountered for event")
        row = str(i)+' N/A N/A N/A N/A N/A N/A N/A N/A N/A N/A N/A \n'
        work = False

    if (work == True): #Traces for event have been succesfully read so lets try to measure splittiing!
        arrival = model.get_travel_times(st[0].stats.sac.evdp,st[0].stats.sac.gcarc,["SKS"]) # Calculates expected SKS arrival for the current event
        time = st[0].stats.starttime # Start time of stream for event. This should be the event origin time.
        print(time) # Prints origin time.
        st.filter("bandpass",freqmin=0.01, freqmax= 0.5,corners=2,zerophase=True)
        SKS = arrival[0].time
        t0 = obspy.core.utcdatetime.UTCDateTime(st[0].stats.starttime + SKS) # SKS arrival time relative to trace start as a UTCDateTime object
        st.trim(0,t0+300)

        east = st[0].data #East compenent seismogram, ready to be put into a Pair object.
        north = st[1].data #North compenent seismogram, ready to be put into a Pair object.
        sample_interval = st[0].stats.delta
        split_pair = sw.Pair(north,east,delta = sample_interval) # Creates the Pair object
        #split_pair.set_window(t0-30,t0+30)
        print("SKS arrival time for event","is ",SKS)
        split_pair.plot(pick=True,marker=SKS) # Plots traces and particle motion, with a marker for the SKS arrival. Plot window also allows manual picking of the window.
        window1 = split_pair.wbeg()
        window2 = split_pair.wend()
        print("Window Begins at",window1,"and ends at",window2)
        split = sw.EigenM(split_pair,lags=(4,) )
        split.save(filename)


        date = str(time.year)+" "+ str(time.month).zfill(2)+" "+str(time.day).zfill(2) +" " # Format date information so that is inteligible when output
        t = str(time.hour).zfill(2)+" "+str(time.minute).zfill(2)+" "+str(time.second).zfill(2) # Formating time infomation for printing
        row = str(i)+' '+date+' '+t+' '+str(stat)+' '+str(split.fast)+' '+str(split.dfast)+' '+str(split.tlag)+' '+str(split.dtlag)+'\n' #Row of data to be written to output textfile
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
