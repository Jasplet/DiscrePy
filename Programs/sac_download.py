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
import timeit
import sys
import numpy as np
import obspy.core.utcdatetime
import os.path
from obspy.clients.fdsn.client import Client
from obspy.clients.fdsn.header import (FDSNException)
from obspy.clients.fdsn.header import (FDSNNoDataException)
from obspy.core import AttribDict
from obspy.clients import iris
import time
from multiprocessing import Pool, current_process
from functools import partial
import contextlib
####################################
def main(mode,outdir,event_list=None,stat_list=None,batch=False):
    """
    Main Routine for this program. Holds highest level logic

    event_list - [str] path to the event list that you want to download data for IF running in batch mode.

    batch - [bool] True/False switch. If true then sac_download will loop over all stations detected in the event list attempts to downlaod data for them

    The event list must be a space delimted text file with the columns for station code, event time and even lat/long in the following format:
    STAT DATE TIME EVLA EVLO
    EX01 1900001 1234 56.78 90.12

    Note that Date has the format yyyyjjj and TIME has the format hhmm. Sac_download will search for the event and add in seconds to origin time
    """
    #print('Have you checked the downloaded streams file names you want to use??? [10 s pause]')
    ext = input('Enter unique identifier for the _downloaded_streams files: ')
    print('Ok now lets get started')
    start = time.time()
    # Reads our event list as a pandas datatframe
    # The converters kwarg fr TIME will stop pandas from stripping off the leading zeros (but time is now a string)

    df = pd.read_csv(event_list,delim_whitespace=True,converters={'TIME': lambda x: str(x)})
    (attempts,dwn,fdsnx,ts,ex) = 0,0,0,0,0
    sum = []
    if mode is 'single':

        stations = df.STAT.unique() # Identify the unique stations provided in the event staton list
        if batch is True:
            # for station in stations:
            dwn = partial(run_download,df=df,ext=ext,out=outdir) # make partial fucntion so we can sue map to parrallise

            with contextlib.closing( Pool(processes = 4) ) as pool:
            #           Iterate over stations in the station list.
                pool.map(dwn, stations) # Map list of unique stations to downloader

                # (a,d,fd,t,x,s)= run_download(df=df,station=station,ext=ext,out=outdir)
                # attempts += a # counter for number of attmepts
                # dwn += d # count number of downloads
                # fdsnx += fd # coutn number of fdsn exceptions
                # ts += t # count numbers of streams that were too short (ideally 0)
                # ex += x # numbers of files that already existed
                # # sum.append(s[0])

        elif batch is False:
            station = input('Input Station Name > ')
            (attempts,dwn,fdsnx,ts,ex,sum) = run_download(df,station,ext,outdir)
    elif mode is 'sep':
        #If event and station lists are seperate
        #We dont need to filter df by station and instead can pass it straight in
        sdf= pd.read_csv(stat_list,delim_whitespace=True)
        # print(sdf)
        stations = sdf.STAT
        for station in stations:

            (a,d,fd,t,x,s)= run_download(df,station,ext,outdir,sep=True,sdf=sdf[sdf['STAT'] == station])
            attempts += a
            dwn += d
            fdsnx += fd
            ts += t
            ex += x
            # print(s)
            sum.append(s[0])
    # print('{:03d} download attempts were made, {:02d} were successful, {:02d} hit FDSNNoDataExceptions, {:02} were incomplete and {:02d} have already been downloaded'.format(attempts,dwn,fdsnx,ts,ex))
    print(sum)
    with open('{}/{}_downloaded_streams.txt'.format(outdir,outdir.split('/')[-1]),'w+') as writer:
        [writer.write('{}\n'.format(id)) for id in sum]

    end = time.time()
    runtime = end - start
    print('Runtime was {}\n'.format(runtime))


def run_download(station,df,ext,out,sep=False,sdf=None):
    """
    Function that runs the downloading process for a given station (or stations)
    """
    data = df[(df['STAT'] == station)]
    data = data.reset_index()
    del data['index']
    (attempts,dwn,fdsnx,ts,ex) = 0,0,0,0,0
    sum = []
    print('Begin download of data from station {}'.format(station))
    Instance = Downloader(data,station,out)
    if Instance.attempts == 0:
        ''' i.e is this the first attempt?  '''
        # print(Instance.attempts)
        Instance.outfile = open('/{}/{}/{}_downloaded_streams_{}.txt'.format(out,station,station,ext),'w+')

    stat_found = Instance.download_station_data()
    if stat_found is True:
        # print(Instance.data)
        for i in range(0,len(Instance.data)):
            # print(station, Instance.data.DATE[i])
        #Loop over events for the given station Instance

            Instance.set_event_data(i,sep)
            print('Station: {} Date: {}, Time: {}'.format(station,Instance.date,Instance.time))
            for channel in ['BHN','BHE','BHZ']:

                Instance.download_traces(channel)

    else:
        print('Station {} could not be found'.format(station))

    Instance.outfile.close()
    print('Dowloads for station {} complete, {} events downloaded'.format(station,Instance.dwn))
    print(Instance.attempts,Instance.dwn,Instance.fdsnx,Instance.ts,Instance.ex,Instance.summary)
    return(Instance.attempts,Instance.dwn,Instance.fdsnx,Instance.ts,Instance.ex,Instance.summary)

class Downloader:

    def __init__(self,df,station,outdir):

        self.station = station
        self.data = df
        self.out = outdir
        self.summary = [] # list to hold all tr_ids
        # print(self.data)
#           Resets indexing of DataFrame

        # print('{}/{}_downloaded_streams.txt'.format(outdir,outdir.split('/')[-1]))
        try:
            #print('Make /Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}'.format(station))
            os.mkdir('{}/{}'.format(self.out,station))
        except FileExistsError:
            print('It already exists, Hooray! Less work for me!')

            # pass
    #   Made

        #self.outfile = open('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}_downloaded_streams_Jacks_Split.txt'.format(station,station),'w+')

        self.attempts = 0 #Counter for how many attempted downloads there were
        self.fdsnx = 0 #Counter for how many attempts hit a FDSNNoDataException
        self.dwn = 0 #Counter for how many events were downloaded
        self.ex = 0 #Counter for how many event already exist in filesystem and therefore werent downloaded
        self.ts = 0 #Counter for events who;s traces are too short.
        self.fdsnclient_evt = Client('IRIS') # Serparate client for events (hopefully to get round the "no event avialbel bug")
        self.fdsnclient = Client('IRIS')
#       Download Station Data

    def download_station_data(self):
        """
        Download or read important station data and make sure it is right
        """
        try:
            stat =  self.fdsnclient.get_stations(channel='BH?',station='{}'.format(self.station))
            self.networks = stat.networks
            self.stla = stat.networks[0].stations[0].latitude
            self.stlo = stat.networks[0].stations[0].longitude
            # print(self.network)
            if len(stat.networks) > 1:
                print('More than one network code found!')
            return True
        except FDSNNoDataException:
            return False

    def set_event_data(self,i,sep):
        """
        Function to download event information so we can get mroe accurate start times
        """
        self.evla = self.data.EVLA[i]
        self.evlo = self.data.EVLO[i]

        self.date = self.data.DATE[i]
        if 'TIME' in self.data.columns:
            self.time = self.data.TIME[i]
        else:
            self.time = '0000'
#       get windows (if they exist). we will need these when generating the sac files
        self.wbeg = self.data.WBEG[i]
        self.wend = self.data.WEND[i]

        datetime = str(self.date) + "T" + self.time #Combined date and time inputs for converstion t UTCDateTime object
        self.start = obspy.core.UTCDateTime(datetime)

        try:
            if 'TIME' in self.data.columns:
                end = self.start + 60
                print('Search starts {} , ends at {}'.format(self.start,end))
                cat = self.fdsnclient_evt.get_events(starttime=self.start,endtime=self.start+86400 ,latitude=self.evla,longitude=self.evlo,maxradius=0.25,minmag=5.5) #Get event in order to get more accurate event times.
                # self.time = '{:02d}{:02d}{:02d}'.format(cat[0].origins[0].time.hour,cat[0].origins[0].time.minute,cat[0].origins[0].time.second)
            else:
                # No Time so we need to search over the whole day
                end = self.start + 86400

                print('Search starts {} , ends at {}'.format(self.start,end))
                cat = self.fdsnclient_evt.get_events(starttime=self.start,endtime=self.start+86400 ,latitude=self.evla,longitude=self.evlo,maxradius=0.25,minmag=5.5) #Get event in order to get more accurate event times.
            if len(cat) > 1:
                print("WARNING: MORE THAN ONE EVENT OCCURS WITHIN 5km Search!!")
                print('Selecting Event with the largest magnitude')
                # Select biggest magnitude
                max_mag = max([cat[j].magnitudes[0].mag for j in [i for i,c in enumerate(cat)]])
                cat = cat.filter('magnitude >= {}'.format(max_mag))
                # print(cat)

            self.date = '{:04d}{:03d}'.format(cat[0].origins[0].time.year,cat[0].origins[0].time.julday)
            self.time = '{:02d}{:02d}{:02d}'.format(cat[0].origins[0].time.hour,cat[0].origins[0].time.minute,cat[0].origins[0].time.second)

            self.start.minute = cat[0].origins[0].time.minute
            self.start.hour = cat[0].origins[0].time.hour
            # print(self.time)

            self.start.second = cat[0].origins[0].time.second

            # Lines commented out as they are only needed if TIME is prvoided as hhmm (For Deng's events there is
            # no TIME provided so we just have to used the event time downloaded)
            # if self.start.minute != cat[0].origins[0].time.minute:
            #     self.time = self.time[:2] + str(cat[0].origins[0].time.minute) # Time is hhmm so we subtract the old minute value and add the new one

            dep = cat[0].origins[0].depth
            if dep is not None:
                self.evdp = dep/1000.0 # divide by 1000 to convert depth to [km[]
            else:
                self.evdp = 10.0 #Hard code depth to 10.0 km if evdp cannot be found
        except FDSNNoDataException:
            print("No Event Data Available")
            self.evdp = 0
        except FDSNException:
            print("FDSNException for get_events")
            # pass
        # elif sep is True:
        #     self.start = obspy.core.UTCDateTime('{}'.format(self.data.DATE[i])) #iso8601=True
        #     self.date = '{:04d}{:03d}'.format(self.start.year,self.start.julday)
        #     self.time = '{:02d}{:02d}{:02d}'.format(self.start.hour,self.start.minute,self.start.second)
        #     self.evdp = self.data.EVDP[i]

    def download_traces(self,ch,n=0):
        """
        Function that downloads the traces for a given event and station
        """
        # if len(self.time) is 6:
        print('Start: {}. self.time: {}'.format(self.start,self.time))
        tr_id = "{}/{}/{}_{}_{}_{}.sac".format(self.out,self.station,self.station,self.date,self.time,ch)
        print('Network code is {}, n is {}'.format(self.networks[n].code,n))
        # elif len(self.time) is 4:
            # tr_id = "{}/{}/{}_{}_{}{}_{}.sac".format(self.out,self.station,self.station,self.date,self.time,self.start.second,ch)
        # print("Looking for :", tr_id)


        if ch == 'BHE':
            self.attempts += 1 # Counts the number of traces that downloads are attempted for

        if os.path.isfile(tr_id) == True:
            print("{} exists. It was not downloaded".format(tr_id)) # File does not exist

            if ch == 'BHE':
                out_id = '_'.join(tr_id.split('_')[0:-1])
                self.outfile.write('{}_\n'.format(out_id))
                self.summary.append(out_id)
                self.ex += 1
        else:
            # print("It doesnt exists. Download attempted")
            st = obspy.core.stream.Stream() # Initialises our stream variable

            if self.networks[n] is 'BK':
                download_client = obspy.clients.fdsn.Client('NCEDC')
            else:
                download_client = obspy.clients.fdsn.Client('IRIS')
            try:
                st = download_client.get_waveforms(self.networks[n].code,self.station,'*',ch,self.start,self.start + 3000,attach_response=True)
                print(len(st))
                if len(st) > 1 :
                    print("WARNING: Unxecpected number of traces {} (more than one trace per channel) downloaded for event tr_id {}, skipping".format(len(st),tr_id))
                    
                else:
                    dist_client = iris.Client() # Creates client to calculate event - station distance
                    print('STLA {} STLO {} EVLA {} EVLO {}'.format(self.stla,self.stlo,self.evla,self.evlo))
                    self.d = dist_client.distaz(stalat=self.stla,stalon=self.stlo,evtlat=self.evla,evtlon=self.evlo)
                    print('Source-Reciever distance is {}'.format(self.d['distance']))
                    # if (self.d['distance'] >= 85.0) or (self.d['distance'] >=145.0): # For SKS, SKKS data
                    if (self.d['distance'] >= 50.0) or (self.d['distance'] >=85.0): # For ScS data
                        if st[0].stats.endtime - st[0].stats.starttime >= 2000:
                            # print('Record length is {}, which is ok'.format(st[0].stats.endtime - st[0].stats.starttime))
                            self.write_st(st,tr_id)
                            if ch == 'BHE':
                                self.dwn += 1
                                out_id = '_'.join(tr_id.split('_')[0:-1])
                                self.outfile.write('{}_\n'.format(out_id))
                                self.summary.append(out_id)

                        else:
                            print('Record length is {}, which is too short'.format(st[0].stats.endtime - st[0].stats.starttime))
                            if ch == 'BHE':
                                self.ts += 1
                    else:
                        print("Source Reciever Distance is too small")
                        if ch == 'BHE':
                            self.ts += 1

            except FDSNException:
                print(n+1, len(self.networks))
                if (len(self.networks) > 1) and (n+1 < len(self.networks)):
                    # Here we are testing if there is more than one network code (becuase if theres only 1 then this is no use
                    # and if there are any untried network codes
                    n +=1 # counter for network ID use
                    # Try to download again
                    print('No data, but multiple network codes. So trying again (for the {}th time)'.format(n))
                    self.download_traces(ch,n)

                print('No Data Exception??')
                if ch == 'BHE':
                    self.fdsnx += 1

    def write_st(self,st,tr_id):
        """

        """
        # print('Writing {}'.format(tr_id))
        st[0].write('holder_{}.sac'.format(current_process().pid), format='SAC',) # Writes traces as SAC files
        #st.plot()
        st_2 = obspy.core.read('holder_{}.sac'.format(current_process().pid))
        #sac = AttribDict() # Creates a dictionary sacd to contain all the header information I want.
        ## Set origin times
        st_2[0].stats.sac.nzyear = self.start.year
        st_2[0].stats.sac.nzjday = self.start.julday
        st_2[0].stats.sac.nzhour = self.start.hour
        st_2[0].stats.sac.nzmin = self.start.minute
        st_2[0].stats.sac.nzsec = self.start.second
        st_2[0].stats.sac.nzmsec = self.start.microsecond
        ## Station Paramters
        st_2[0].stats.sac.stla = self.stla
        st_2[0].stats.sac.stlo = self.stlo
        ## Event Paramters
        st_2[0].stats.sac.evla = self.evla#cat[0].origins[0].latitude # Event latitude
        st_2[0].stats.sac.evlo = self.evlo#cat[0].origins[0].longitude # Event longitude
        st_2[0].stats.sac.evdp = self.evdp#cat[0].origins[0].depth/1000 # Event depth
        st_2[0].stats.sac.kstnm = '{:>8}'.format(self.station)
        # print('stla = {}, stlo = {}, evla = {}, evlo = {}'.format(stla,stlo,evla,evlo))


        st_2[0].stats.sac.gcarc = self.d['distance'] # d.values returns the values from dictionary d produced by distaz. list converts this to a list attribute which can then be indexed to extract the great cricle distance in degrees
        st_2[0].stats.sac.dist = self.d['distancemeters']/1000 # Distnace in kilometers
        st_2[0].stats.sac.baz = self.d['backazimuth'] # Backzimuth (Reciever - SOurce)
        st_2[0].stats.sac.az = self.d['azimuth'] # Azimuth (Source - Receiver)
        st_2[0].write(tr_id, format='SAC',byteorder=1)

        # Perturb the pre-existing windows slightly (as they are for long-period data)
        # This is so we can use nwind in SHEBA later...
        st_2[0].stats.sac.user0 = self.wbeg - 2.5
        st_2[0].stats.sac.user1 = self.wbeg + 2.5
        st_2[0].stats.sac.user2 = self.wbeg - 2.5
        st_2[0].stats.sac.user3 = self.wend + 2.5

if __name__ == '__main__':
    # This block allows this mess of code to be run outside of ipython as a script
    start = timeit.default_timer()
    print('You are running sac_download from the command line. I am assuming there is a signle event list (e.g. RDB, SDB etc) list to work from')
    print('If you have seperate event and station lists, a) Why? b) You will have to use ipython')
    out = sys.argv[1]
    evts = sys.argv[2]

    main(mode='single',outdir=out,event_list=evts,batch=True)

    ######################################################################################
    # End Timing of run
    ######################################################################################
    end = timeit.default_timer()
    runtime = end - start
    print('The runtime of main is {} minutes'.format(runtime/60))
