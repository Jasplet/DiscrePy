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
def trace_download(dates,times,evla,evlo,evdp,stla,stlo,station,network,type, n_events):
    ## Function to download and save traces for a pre-determined set of events

    too_short = 0
    no_dat = 0
    outfile = open('/Users/ja17375/Python/SKS_Splitting/Data/{}/{}_downloaded_{}_streams.txt'.format(station,station,type),'w+')
    for i in range(0, n_events):
        # Test if trace is already downloaded and create a UTCDateTime object for the Starttime
        ps = str(dates[i]) + "T" + str(times[i]).zfill(4)
        start = obspy.core.UTCDateTime(ps) #iso8601=True
        # client = obspy.clients.fdsn.Client("IRIS")
        n = 0
        client = obspy.clients.fdsn.Client('IRIS') # BKS Station data is hosted by NCEDC so we need this as the client
        try:
            cat = client.get_events(starttime=start-60,endtime=start+60 ,latitude=evla[i],longitude=evlo[i],maxradius=0.5) #Get event in order to get more accurate event times.
            if len(cat) > 1:
                print("WARNING: MORE THAN ONE EVENT OCCURS WITHIN 5km Search!!")

            start.microsecond = cat[0].origins[0].time.microsecond
            start.second = cat[0].origins[0].time.second
            start.minute = cat[0].origins[0].time.minute
            sec = str(start.second).zfill(2)
            n = 1
        except FDSNNoDataException:
            print("No Event Data Available")
        except FDSNException:
            print("FDSNException for get_events")

        channel = ["BHN","BHZ","BHE"]

        date = '{:07d}'.format(dates[i])
        if n == 0:
            time = '{:04d}00'.format(times[i])
        elif n ==1:
            time = '{:04d}{}'.format(times[i],sec)

        for ch in channel:

            id_tst = "/Users/ja17375/Python/SKS_Splitting/Data/{}/{}_{}_{}_{}.sac".format(station,station,date,time,ch)
            print("Looking for :", id_tst)
            if os.path.isfile(id_tst) == True:
                print("It exists. It was not downloaded") # File does not exist
            else:
                print("It does not exist. Download attempted")
                st = obspy.core.stream.Stream() # Initialises our stream variable

                loc = '--'

                #print("Event: ", i, ". Station: ", station, ". Channel: ", ch)
                # test_00 = obspy.core.UTCDateTime("2011-05-02T19:07:00") # This Criterion is ONLY for AAM
                # if start >= test_00:
                #  # Location "00" for AAM went live at 2011/05/02 (122) 19:07:00
                #  # Location "10" is for NEW
                #    loc = "10"
                #
                # elif start < obspy.core.UTCDateTime("2000-01-01T00:00:00"):
                #     loc = "--"
                #
                # else:
                #     loc = "--"

                try:

                    download_client = obspy.clients.fdsn.Client('IRIS')
                    st = download_client.get_waveforms(network,station,'??',ch,start,start + 3000,attach_response=True)

                    if len(st) > 3:
                        print("WARNING: More than three traces downloaded for event ", i)
                    if ((st[0].stats.endtime - st[0].stats.starttime) >= 2999.0):

                        tr_id = "/Users/ja17375/Python/SKS_Splitting/Data/{}/{}_{}_{}_{}.sac".format(station,station,date,time,ch)
                        st[0].write('holder.sac', format='SAC',) # Writes traces as SAC files
                        #st.plot()
                        st_2 = obspy.core.read('holder.sac')
                        #sac = AttribDict() # Creates a dictionary sacd to contain all the header information I want.
                        ## Station Paramters
                        st_2[0].stats.sac.stla = stla
                        st_2[0].stats.sac.stlo = stlo
                        ## Event Paramters
                        st_2[0].stats.sac.evla = evla[i]#cat[0].origins[0].latitude # Event latitude
                        st_2[0].stats.sac.evlo = evlo[i]#cat[0].origins[0].longitude # Event longitude
                        st_2[0].stats.sac.evdp = evdp[i]#cat[0].origins[0].depth/1000 # Event depth
                        dist_client = iris.client.Client() # Creates client to calculate event - station distance
                        d = dist_client.distaz(stalat=stla,stalon=stlo,evtlat=evla[i],evtlon=evlo[i])
                        st_2[0].stats.sac.gcarc = d['distance'] # d.values returns the values from dictionary d produced by distaz. list converts this to a list attribute which can then be indexed to extract the great cricle distance in degrees
                        st_2[0].stats.sac.dist = d['distancemeters']/1000 # Distnace in kilometers
                        st_2[0].stats.sac.baz = d['backazimuth'] # Backzimuth (Reciever - SOurce)
                        st_2[0].stats.sac.az = d['azimuth'] # Azimuth (Source - Receiver)
                        ## File Type
                        #sacd.iftype = 1
                        #print(st_2[0].stats)
                        #st[0].stats.sac = sac
                        #st_2.plot()
                        st_2[0].write(tr_id, format='SAC')
                        print("The trace ", tr_id, "was downloaded and saved!")

                        if ch == 'BHE':
                            outfile.write('{}\n'.format(tr_id[0:-7]))


                    else:
                        print("Trace is too short.")
                        too_short = too_short+1
                except FDSNException:
                    #print("No data available for event", i+1, ":-(")
                    no_dat = no_dat + 1 #Counter for number of events where no data was availables

def split_read(station,network,autoqc):
    """
    Initialises some variable and call the trace_dowload function for a given station
    """
    data = pd.read_csv('/Users/ja17375/Python/SKS_Splitting/Data/Jacks_SKS_RAW.txt',delim_whitespace=True,keep_default_na=False)
    num_splits = data[data.AUTOQC == 'split' ]['STAT'].value_counts()
    # stat=num_splits.index[0] ## When made into a function this will be the requested variable (or one of them at least)

    split1 = data[(data['STAT'] == station) & (data['AUTOQC'] == autoqc)]
    split1 = split1.reset_index()
    del split1['index']
    ### This appear to be the best way to parse the data that we want
    # Resets the indexing of the rows to make them easier to find


    date = split1.DATE
    time = split1.TIME
    dist = split1.DIST
    evdp = split1.EVDP
    evla = split1.EVLA
    evlo = split1.EVLO
    stla = split1.STLA[0]
    stlo = split1.STLO[0]

    no_events = len(date) # Each event has a date, so the length of date gives the number of events

    # for i in range(0 , split1.shape[0]): # split1.shape[0] returns the number of rows in split1
    # Call the function trace_download, to download traces
    trace_download(date,time,evla,evlo,evdp,stla,stlo,station,network,autoqc, no_events)
