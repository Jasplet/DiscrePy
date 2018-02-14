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
def split_read(station,network='*'):
    """
    Initialises some variable and call the trace_dowload function for a given station
    """
    raw = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Python/Data/Jacks_SKS_RAW.txt',delim_whitespace=True,keep_default_na=False)

    data = raw[(raw['STAT'] == station) & (raw['AUTOQC'] != 'fail')]
    ### This appear to be the best way to parse the data that we want
    # Resets the indexing of the rows to make them easier to find
    data = data.reset_index()
    del data['index']

#   Make target directory. the exist_ok=True flag means that if it already exists then there will be no error
    try:
        print('Make /Users/ja17375/Shear_Wave_Splitting/Python/Data/{}'.format(station))
        os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Python/Data/{}'.format(station))
    except FileExistsError:
        print('It already exists, Hooray! Less work for me!')
#   Made
    outfile = open('/Users/ja17375/Shear_Wave_Splitting/Python/Data/{}/{}_downloaded_streams.txt'.format(station,station),'w+')


    attempts = 0 #Counter for how many attempted downloads there were
    fdsnx = 0 #Counter for how many attempts hit a FDSNNoDataException
    dwn = 0 #Counter for how many events were downloaded
    ex = 0 #Counter for how many event already exist in filesystem and therefore werent downloaded
    ts = 0

    for i in range(0,len(data)):
    # Call the function trace_download, to download traces
        dwn,fdsnx, ex, ts = trace_download(data.DATE[i],data.TIME[i],data.EVLA[i],data.EVLO[i],data.EVDP[i],data.STLA[i],data.STLO[i],station,network,outfile,fdsnx,ex,dwn,ts)
        attempts += 1

    print('{:03d} download attempts were made, {:02d} were successful, {:02d} hit FDSNNoDataExceptions, {:02} were incomplete and {:02d} aboarted as the data has already been downloaded'.format(attempts,dwn,fdsnx,ts,ex))

def trace_download(date,time,evla,evlo,evdp,stla,stlo,station,network,outfile,fdsnx,ex,dwn,ts):
    ## Function to download and save traces for a pre-determined set of events


    datetime = str(date) + "T" + str(time).zfill(4) #Combined date and time inputs for converstion t UTCDateTime object
    start = obspy.core.UTCDateTime(datetime) #iso8601=True

    client = obspy.clients.fdsn.Client('IRIS') #
    try:
        cat = client.get_events(starttime=start-60,endtime=start+60 ,latitude=evla,longitude=evlo,maxradius=0.5) #Get event in order to get more accurate event times.
        if len(cat) > 1:
            print("WARNING: MORE THAN ONE EVENT OCCURS WITHIN 5km Search!!")

        start.second = cat[0].origins[0].time.second

        if start.minute != cat[0].origins[0].time.minute:
            time = (time - start.minute) + cat[0].origins[0].time.minute # Time is hhmm so we subtract the old minute value and add the new one

    except FDSNNoDataException:
        print("No Event Data Available")
    except FDSNException:
        print("FDSNException for get_events")


    channel = ["BHN","BHZ","BHE"]
    for ch in channel:

        tr_id = "/Users/ja17375/Shear_Wave_Splitting/Python/Data/{}/{}_{:07d}_{:04d}{:02d}_{}.sac".format(station,station,date,time,start.second,ch)
        # print("Looking for :", id_tst)
        if os.path.isfile(tr_id) == True:
            print("It exists. It was not downloaded") # File does not exist
            if ch == 'BHE':
                outfile.write('{}\n'.format(tr_id[0:-7]))
                ex += 1
        else:
            # print("It does not exist. Download attempted")
            st = obspy.core.stream.Stream() # Initialises our stream variable
            try:
                if network is 'BK':
                    download_client = obspy.clients.fdsn.Client('NCEDC')
                else:
                    download_client = obspy.clients.fdsn.Client('IRIS')

                st = download_client.get_waveforms(network,station,'??',ch,start,start + 3000,attach_response=True)

                if len(st) > 3:
                    print("WARNING: More than three traces downloaded for event ", tr_id)
                if ((st[0].stats.endtime - st[0].stats.starttime) >= 2999.0):


                    st[0].write('holder.sac', format='SAC',) # Writes traces as SAC files
                    #st.plot()
                    st_2 = obspy.core.read('holder.sac')
                    #sac = AttribDict() # Creates a dictionary sacd to contain all the header information I want.
                    ## Station Paramters
                    st_2[0].stats.sac.stla = stla
                    st_2[0].stats.sac.stlo = stlo
                    ## Event Paramters
                    st_2[0].stats.sac.evla = evla#cat[0].origins[0].latitude # Event latitude
                    st_2[0].stats.sac.evlo = evlo#cat[0].origins[0].longitude # Event longitude
                    st_2[0].stats.sac.evdp = evdp#cat[0].origins[0].depth/1000 # Event depth
                    dist_client = iris.Client() # Creates client to calculate event - station distance
                    # print('stla = {}, stlo = {}, evla = {}, evlo = {}'.format(stla,stlo,evla,evlo))

                    d = dist_client.distaz(stalat=stla,stalon=stlo,evtlat=evla,evtlon=evlo)

                    st_2[0].stats.sac.gcarc = d['distance'] # d.values returns the values from dictionary d produced by distaz. list converts this to a list attribute which can then be indexed to extract the great cricle distance in degrees
                    st_2[0].stats.sac.dist = d['distancemeters']/1000 # Distnace in kilometers
                    st_2[0].stats.sac.baz = d['backazimuth'] # Backzimuth (Reciever - SOurce)
                    st_2[0].stats.sac.az = d['azimuth'] # Azimuth (Source - Receiver)


                    ## File Type
                    #sacd.iftype = 1
                    #print(st_2[0].stats)
                    #st[0].stats.sac = sac
                    #st_2.plot()
                    st_2[0].write(tr_id, format='SAC',byteorder=1)
                    # print("The trace ", tr_id, "was downloaded and saved!")
                    dwn += 1
                    if ch == 'BHE':
                        outfile.write('{}\n'.format(tr_id[0:-7]))


                else:
                    print("Trace is too short.")
                    ts += 1
            except FDSNException:
                if ch == 'BHE':
                    fdsnx += 1

    return dwn, fdsnx, ex,ts
