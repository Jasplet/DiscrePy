#! /usr/bin/env python
##############################
#   Program: interact_sheba.py
#
##############################
#   Author: J. Asplet
##############################
#   A python program to act as a wrapper to sheba.
#   It read, filters, and trims SAC files from a
#   provided list of events and passes each event
#   to sheba for shear wave splitting analysis.
#   When the analysis is complete the program will
#   collect the output files for each event and it
#   will concatenate them together into a summary
#   output file!
#
##############################
#   Import Statements
##############################
#   Standard Packages - all freely available
import obspy as ob
import numpy as np
import pandas as pd
import matplotlib.pyplot as ply #Just incase
import subprocess as sub
import os.path
##############################
#   Import other scripts in Programs/
import Split_Read as sr
import Split_Measure as sm #Just in incase
import plot as sp
##############################
def main(phase='SKS'):
    """
    Main - call this function to run the interface to sheba
    """
   #First Indentify the possible station that we could have data for
   #This way we know what directory paths to look in for the sac files
    path = '/Users/ja17375/Shear_Wave_Splitting'
    statlist ='{}/Data/StationList.txt'.format(path)
    stations = pd.read_csv(statlist,delim_whitespace=True).STAT
    i = 1 # Counter
   #Loop over all stations in the list.
    for station in stations:
        #Each station SHOULD have its own directory within Data/SAC_files
        #If the data has been downloaded. So lets look for directorys that exist
        dir_path = '{}/Data/SAC_files/{}'.format(path,station)
        if os.path.isdir(dir_path):
            #Happy Days! The data directory exists!
            with open('{}/{}_downloaded_streams.txt'.format(dir_path,station),'r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
                for line in reader.readlines():
                    line.strip('\n')
                    st_id = '{}BH?.sac'.format(str(line[0:-1]))
                    st = ob.read(st_id)
                    if len(st) is 3:
                        Event = Interface(st)
                        if phase is 'SKS':
                            Event.process(station,phase,i=i,path='{}/Sheba/SAC'.format(path))
                            Event.sheba(station,phase,i=i,path='{}/Sheba/SAC'.format(path))
                            i +=1
                        elif phase is 'SKKS':
#                           When we look for SKKS the traces will need to be trimmmed at different times!
                            Event.process(t1 = 1800, t2 = 2000)
                            i +=1
                    else:
                        pass
        else:
            print('The directory {}/Data/SAC_files/{} does not exists'.format(path,station))



class Interface:
    """
    Class which will act as the interface to sheba.
    The "subprocess" sheba will be a bound method
    """
    def __init__(self,st):
        # self.date = date
        # self.time = time
        self.BHE = st.select(channel='BHE')
        self.BHE[0].stats.sac.cmpinc = 90
        self.BHE[0].stats.sac.cmpaz = 90
        self.BHN = st.select(channel='BHN')
        self.BHN[0].stats.sac.cmpinc = 90
        self.BHN[0].stats.sac.cmpaz = 0
        self.BHZ = st.select(channel='BHZ')
        self.BHZ[0].stats.sac.cmpinc = 0
        self.BHZ[0].stats.sac.cmpaz = 0

    def process(self,station,phase,c1=0.01,c2=0.5,t1=1400,t2=1600,i=None,path=None):
        """
        Function to bandpass filter and trim the components
        t1 - [s] Lower bound of trim, time relative to event time
        t2 - [s] Upper bound of trim, time relative to event time
        c1 - [Hz] Lower corner frequency
        c2 - [Hz] Upper corner frequency
        By default traces will trim between 1400 - 1600s and will be filtered between 0.01Hz-0.5Hz
        """


#       De-mean and detrend each component
        self.BHN.detrend(type='demean') #demeans the component
        self.BHE.detrend(type='demean')
        self.BHZ.detrend(type='demean')
#       Detrend
        self.BHN.detrend(type='simple') #De-trends component
        self.BHE.detrend(type='simple') #De-trends component
        self.BHZ.detrend(type='simple') #De-trends component
#       Filter each component. Bandpass flag gives a bandpass-butterworth filter
        self.BHN.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHE.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHZ.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
#       Now trim each component to the input length (default is 1400,1600)
        self.BHN.trim(self.BHN[0].stats.starttime + t1,self.BHN[0].stats.starttime + t2)
        self.BHE.trim(self.BHE[0].stats.starttime + t1,self.BHE[0].stats.starttime + t2)
        self.BHZ.trim(self.BHZ[0].stats.starttime + t1,self.BHZ[0].stats.starttime + t2)
#       Add windowing ranges to sac headers user0,user1,user2,user3 [start1,start2,end1,end2]
#       As we will use Teanby's multiwindowing subroutine, my window ranges will be at the start, the middle and the end of the trim (1400,1500,1500,1600 for SKS)
        self.BHN[0].stats.sac.user0,self.BHN[0].stats.sac.user1,self.BHN[0].stats.sac.user2,self.BHN[0].stats.sac.user3 = (t1,(t1+t2/2),(t1+t2/2),t2)
        self.BHE[0].stats.sac.user0,self.BHE[0].stats.sac.user1,self.BHE[0].stats.sac.user2,self.BHE[0].stats.sac.user3 = (t1,(t1+t2/2),(t1+t2/2),t2)
        self.BHZ[0].stats.sac.user0,self.BHZ[0].stats.sac.user1,self.BHZ[0].stats.sac.user2,self.BHZ[0].stats.sac.user3 = (t1,(t1+t2/2),(t1+t2/2),t2)
#       Now write out the three processed components
#       Naming depends on whether this is being executed as a test or within a loop
#       where a counter should be provided to prevent overwriting.
        if i is not None:
            self.BHN.write('{}/{}_{}_{}.BHN'.format(path,station,phase,i),format='SAC',byteorder=1)
            self.BHE.write('{}/{}_{}_{}.BHE'.format(path,station,phase,i),format='SAC',byteorder=1)
            self.BHZ.write('{}/{}_{}_{}.BHZ'.format(path,station,phase,i),format='SAC',byteorder=1)
        else:
            self.BHN.write('{}.BHN'.format(station),format='SAC',byteorder=1)
            self.BHE.write('{}.BHE'.format(station),format='SAC',byteorder=1)
            self.BHZ.write('{}.BHZ'.format(station),format='SAC',byteorder=1)


    def plot_comp(self):
        """
        Quick Function to plot component together on one seismogram
        """
        st = self.BHN + self.BHE + self.BHZ
        st.plot(type='relative')

    def sheba(self,station,phase,i = None,path=None):
        """
        The big one! This function uses the subprocess module to host sac and then runs sheba as a SAC macro
        """
        p = sub.Popen(['sac'],
                     stdout = sub.PIPE,
                     stdin  = sub.PIPE,
                     stderr = sub.STDOUT,
                     encoding='utf8')

        if i is not None:
            s = '''
            echo on\n
            m sheba file {}/{}_{}_{} plot no nwind 10 10 batch yes
            '''.format(path,station,phase,i)
            print(s)
            out =p.communicate(s)
            # print(out[0])
        else:
            s = '''
            echo on\n
            m sheba file {} plot no nwind 10 10
            '''.format(station)

            out =p.communicate(s)
            print(out[0])
## Psuedo code plan for script
# Read Station list

# for stat in list:
    #request sac file list from data/SAC_files/stat/

    #for event in eventlist:
        #st = ob.read(event*

        #filter stream

        #trim stream to reasonable range (1400- 1600)

        #byte-swap stream
        #set headers for inclination/azimuth (call sac subprocess to do this)

        # sac m sheba
    #end loop over eventlist
#end loop over stationlist

#concatenate sheba result files together - possibly add some more metadata as well

#subroutine for stacking error surfaces
