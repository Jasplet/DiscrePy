#! /usr/bin/env python
## This is a script to read the SKS splitting data from Walpole et al 2014 and manipulate that data.

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
from obspy.clients.fdsn.header import (FDSNException)
from obspy.clients.fdsn.header import (FDSNNoDataException)
from obspy.core import AttribDict
from obspy.clients import iris


##################### Start of Function to Request and Save Traces ##################
def trace_download(dates,times,evla,evlo,stla,stlo,station, n_events):
    ## Function to download and save traces for a pre-determined set of events
    too_short = 0
    no_dat = 0
    for i in range(0 n_events):
        # Test if trace is already downloaded and create a UTCDateTime object for the Starttime
        ps = str(dates[i]) + "T" + str(times[i]).zfill(4)
        start = obspy.core.UTCDateTime(ps) #iso8601=True
        yr = str(start.year)
        mn_tst = str(start.month).zfill(2) # Converts Month to string. zfill ensures it is 2 characters by adding a leading zero if required
        d_tst = str(start.day).zfill(2)   # Day of event, this is to contruct expected
        h_tst = str(start.hour).zfill(2) # Converts Hour to string. zfill ensures it is 2 characters by adding a leading zero if required
        min_tst = str(start.minute).zfill(2)

        client = obspy.clients.fdsn.Client("IRIS")
        try:
            cat = client.get_events(starttime=start-60,endtime=start+60 ,latitude=evla[i],longitude=evlo[i],maxradius=0.5) #Get event in order to get more accurate event times.
            print("i tried. there are this many events", len(cat))
            if len(cat) > 1:
                print("WARNING: MORE THAN ONE EVENT OCCURS WITHIN 5km Search!!")

            start.microsecond = cat[0].origins[0].time.microsecond
            start.second = cat[0].origins[0].time.second
            start.minute = cat[0].origins[0].time.minute
        except FDSNNoDataException:
            print("No Event Data Available")
        except FDSNException:
            print("FDSNException for get_events")
        channel = ["BHN","BHZ","BHE"]
        for ch in channel:
            id_tst = "NEW_"  +  str(i) + "_" + ch + ".sac"
            print("Looking for :", id_tst)
            if os.path.isfile(id_tst) == True:
                print("It exists. It was not downloaded") # File does not exist
            else:
                print("It does not exist. Download attempted")
                st = obspy.core.stream.Stream() # Initialises our stream variable
                client = Client("IRIS")
                #print("Event: ", i, ". Station: ", station, ". Channel: ", ch)
                test_00 = obspy.core.UTCDateTime("2011-05-02T19:07:00") # This Criterion is ONLY for AAM
                if start >= test_00:
                 # Location "00" for AAM went live at 2011/05/02 (122) 19:07:00
                 # Location "10" is for NEW
                   loc = "10"
                   # channel = "BH?"
                elif start < obspy.core.UTCDateTime("2000-01-01T00:00:00"):
                    loc = "--"
                    # channel = "BH?"
                else:
                    loc = "--"
                    # channel = "BH?"
                try:
                    st = client.get_waveforms("US",station,loc,ch,start,start + 3000,attach_response=True)
                    if len(st) > 3:
                        print("WARNING: More than three traces downloaded for event ", i)
                    if ((st[0].stats.endtime - st[0].stats.starttime) >= 2999.0):
                        t = obspy.core.UTCDateTime(st[0].stats.starttime, precision = 0)
                        yr = str(t.year)
                        mon = str(t.month).zfill(2) # Converts Month to string. zfill ensures it is 2 characters by adding a leading zero if required
                        day = str(t.day).zfill(2)   # Day of event, this is to contruct expected
                        hour = str(t.hour).zfill(2) # Converts Hour to string. zfill ensures it is 2 characters by adding a leading zero if required
                        minute = str(t.minute).zfill(2)
                        tr_id = "./Data/" + stat + "_" + str(i) +"_" + ch + ".sac"
                        st[0].write('holder.sac', format='SAC') # Writes traces as SAC files
                        #st.plot()
                        st_2 = obspy.core.read('holder.sac')
                        #sac = AttribDict() # Creates a dictionary sacd to contain all the header information I want.
                        ## Station Paramters
                        st_2[0].stats.sac.stla = stla
                        st_2[0].stats.sac.stlo = stlo
                        ## Event Paramters
                        st_2[0].stats.sac.evla = cat[0].origins[0].latitude # Event latitude
                        st_2[0].stats.sac.evlo = cat[0].origins[0].longitude # Event longitude
                        st_2[0].stats.sac.evdp = cat[0].origins[0].depth/1000 # Event depth
                        dist_client = iris.client.Client() # Creates client to calculate event - station distance
                        d = dist_client.distaz(stalat=stla,stalon=stlo,evtlat=evla[i],evtlon=evlo[i])
                        st_2[0].stats.sac.gcarc = list(d.values())[1] # d.values returns the values from dictionary d produced by distaz. list converts this to a list attribute which can then be indexed to extract the great cricle distance in degrees
                        st_2[0].stats.sac.dist = list(d.values())[2]/1000 # Distnace in kilometers
                        ## File Type
                        #sacd.iftype = 1
                        #print(st_2[0].stats)
                        #st[0].stats.sac = sac
                        #st_2.plot()
                        st_2[0].write(tr_id, format='SAC')
                        print("The trace ", tr_id, "was downloaded and saved!")
                    else:
                        #print("Trace is not complete so it was not saved.")
                        too_short = too_short+1
                except FDSNException:
                    #print("No data available for event", i+1, ":-(")
                    no_dat = no_dat + 1 #Counter for number of events where no data was availables
############## End of Function to Request and Save Traces ####################
############## Start of Main script ##########################################
data = pd.read_csv("./Data/Jacks_SKS_RAW.txt",delim_whitespace=True)
num_splits = data[data.AUTOQC == 'split' ]['STAT'].value_counts()
stat=num_splits.index[0] ## When made into a function this will be the requested variable (or one of them at least)
split1 = data[(data['STAT'] == stat) & (data['AUTOQC'] =="split") ]

### This appear to be the best way to parse the data that we want
split1 = split1.reset_index()
del split1['index'] # Resets the indexing of the rows to make them easier to find
no_events = split1.shape[0]
date = split1.DATE
time = split1.TIME
dist = split1.DIST
depth = split1.EVDP
evla = split1.EVLA
evlo = split1.EVLO
stla = split1.STLA[0]
stlo = split1.STLO[0]

# for i in range(0 , split1.shape[0]): # split1.shape[0] returns the number of rows in split1
# Call the function trace_download, to download traces
trace_download(date,time,evla,evlo,stla,stlo,stat, no_events)
