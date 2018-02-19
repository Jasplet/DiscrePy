#! /usr/bin/env python
####################################
#   Program: sac_download.py
####################################
#   Author: Joseph Asplet
####################################
#   This is a program to download seismograms for a given set
#   of events in the SAC format.
####################################
#   Standard Import Statements
####################################
import obspy as ob
import pandas as pd
import numpy as np
import obspy.core.utcdatetime
import os.path
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import (FDSNException)
from obspy.clients.fdsn.header import (FDSNNoDataException)
from obspy.core import AttribDict
from obspy.clients import iris
####################################
def main(event_list=None,batch=False):
    """
    Main Routine for this program. Holds highest level logic
    """
    if batch is True:
        # Reads our event list as a pandas datatframe
        df = pd.read_csv(eventlist,delim_whitespace=True,converters={'TIME': lambda x: str(x)})
        # The converters kwarg fr TIME will stop pandas from stripping off the leading zeros (but time is now a string)
        for station in df.STAT.unique():
#        for each unique station Code
            Instance = Downloader(df,station)

            for i in range(0,len(Instance.data)):
                #Loop over events for the given station Instance
                Instance.download_event_data()
                for channel in ['BHN','BHE','BHZ']:
                    download_traces(channel)
        print('{:03d} download attempts were made, {:02d} were successful, {:02d} hit FDSNNoDataExceptions, {:02} were incomplete and {:02d} aboarted as the data has already been downloaded'.format(Instance.attempts,Instance.dwn,Instance.fdsnx,Instance.ts,Instance.ex))

class Downloader:

    def __init__(self,df,station):

        self.station = station
        self.data = df[(df['STAT'] == station)]
        self.data = self.data.reset_index()
        del self.data['index']
#           Resets indexing of DataFrame

        try:
            print('Make /Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}'.format(station))
            os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}'.format(station))
        except FileExistsError:
            print('It already exists, Hooray! Less work for me!')
    #   Made
        self.outfile = open('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}_downloaded_streams.txt'.format(station,station),'w+')

        self.attempts = 0 #Counter for how many attempted downloads there were
        self.fdsnx = 0 #Counter for how many attempts hit a FDSNNoDataException
        self.dwn = 0 #Counter for how many events were downloaded
        self.ex = 0 #Counter for how many event already exist in filesystem and therefore werent downloaded
        self.ts = 0 #Counter for events who;s traces are too short.
        self.fdsnclient = Client('IRIS')
#       Download Station Data
        stat =  self.fdsnclient.get_stations(channel='BH?',station='{}'.format(self.station))
        self.network = stat.networks[0].code
        self.stla = stat.networks[0].stations[0].latitude
        self.stlo = stat.netowrks[0].stations[0].longitdue

    def download_station_data(self):
        """
        Download or read important station data and make sure it is right
        """
        stat =  self.fdsnclient.get_stations(channel='BH?',station='{}'.format(self.station))
        self.network = stat.networks[0].code
        self.stla = stat.networks[0].stations[0].latitude
        self.stlo = stat.netowrks[0].stations[0].longitdue

    def download_event_data(self):
        """
        Function to download event information so we can get mroe accurate start times
        """
        datetime = str(self.data.DATE[i]) + "T" + self.data.TIME[i] #Combined date and time inputs for converstion t UTCDateTime object
        self.start = obspy.core.UTCDateTime(datetime) #iso8601=True

        try:
            cat = fdsnclient.get_events(starttime=self.start-60,endtime=self.start+60 ,latitude=self.data.EVLA[i],longitude=self.EVLO[i],maxradius=0.5) #Get event in order to get more accurate event times.
            if len(cat) > 1:
                print("WARNING: MORE THAN ONE EVENT OCCURS WITHIN 5km Search!!")

            self.start.second = cat[0].origins[0].time.second

            if self.start.minute != cat[0].origins[0].time.minute:
                self.time = self.time[:2] + str(cat[0].origins[0].time.minute) # Time is hhmm so we subtract the old minute value and add the new one

        except FDSNNoDataException:
            print("No Event Data Available")
        except FDSNException:
            print("FDSNException for get_events")

    def download_traces(self,ch):
        """

        """
        tr_id = "/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}_{:07d}_{:04d}{:02d}_{}.sac".format(self.station,self.station,self.data.DATE[i],self.data.TIME[i],self.start.second,ch)
        #
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

                    write_st(st)
                    self.dwn += 1
                    if ch == 'BHE':
                        outfile.write('{}\n'.format(tr_id[0:-7]))

                else:
                    print("Trace is too short.")
                    self.ts += 1
            except FDSNException:
                if ch == 'BHE':
                    self.fdsnx += 1

    def write_st(self,st):
        """

        """
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
        st_2[0].stats.sac.kstnm = '{:>8}'.format(station)
        dist_client = iris.Client() # Creates client to calculate event - station distance
        # print('stla = {}, stlo = {}, evla = {}, evlo = {}'.format(stla,stlo,evla,evlo))

        d = dist_client.distaz(stalat=stla,stalon=stlo,evtlat=evla,evtlon=evlo)

        st_2[0].stats.sac.gcarc = d['distance'] # d.values returns the values from dictionary d produced by distaz. list converts this to a list attribute which can then be indexed to extract the great cricle distance in degrees
        st_2[0].stats.sac.dist = d['distancemeters']/1000 # Distnace in kilometers
        st_2[0].stats.sac.baz = d['backazimuth'] # Backzimuth (Reciever - SOurce)
        st_2[0].stats.sac.az = d['azimuth'] # Azimuth (Source - Receiver)
        st_2[0].write(tr_id, format='SAC',byteorder=1)

# PsuedoCode - For developoment

# read list of events as daatframe
#    for unique stations:
#       get network/location code information
#       for events in df.stat=stat:
#           get accurate event time
#               for BHN,BHE.BHZ:
#                   request download of trace of set length
#                   save to holder
#                   read holder and add desired sac headers
#
