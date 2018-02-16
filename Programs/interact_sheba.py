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
from subprocess import CalledProcessError
import os.path
import time
##############################
#   Import other scripts in Programs/
import Split_Read as sr
import Split_Measure as sm #Just in incase
import plot as sp
##############################
def main(phase='SKS',batch=False):
    """
    Main - call this function to run the interface to sheba
    This function depends on the existence of a station list file specified by statlist and you have sac data alreayd downloaded
    Path should point to the directory that you want the sheba processing directories to be stored under

    """
   #First Indentify the possible station that we could have data for
   #This way we know what directory paths to look in for the sac files
    path = '/Users/ja17375/Shear_Wave_Splitting'

   #Loop over all stations in the list.
   #################### Time run #################
    start = time.time()
    ################## Start Process #############
    if batch is True:
#       If processing of data for multiple stations is desired
        statlist ='{}/Data/StationList.txt'.format(path)
        stations = pd.read_csv(statlist,delim_whitespace=True).STAT

        for station in stations:
#           Iterate over stations in the station list.
            run_sheba(path,station,phase)

    elif batch is False:
        station = input('Input Station Name > ')
        run_sheba(path,station,phase)
    end = time.time()
    runtime = end - start
    print('The runtime of main is {} seconds'.format(runtime))

def run_sheba(path,station,phase):
    """
    Function that holds the guts of the workflow for preparing SAC files and running sheba
    """

    #Each station SHOULD have its own directory within Data/SAC_files
    #If the data has been downloaded. So lets look for directorys that exist
    dir_path = '{}/Data/SAC_files/{}'.format(path,station)
    if os.path.isdir(dir_path):
        #Happy Days! The data directory exists!
        i = 1 # Counter
        with open('{}/{}_downloaded_streams.txt'.format(dir_path,station),'r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
            for line in reader.readlines():
                line.strip('\n')
                st_id = '{}BH?.sac'.format(str(line[0:-1]))
                st = ob.read(st_id)
                if len(st) is 3:
                    Event = Interface(st)
                    Event.process(station,phase)
                    outdir = '{}/Sheba/SAC/{}/{}'.format(path,station,phase)

                    try:

                        Event.write_out(station,phase,i=i,path=outdir)
                    except OSError:
                        print('Directory for writing outputs do not all exist. Initialising')
                        os.makedirs(outdir)
                        Event.write_out(station,phase,i=i,path=outdir)

                    Event.sheba(station,phase,i=i,path='{}/Sheba/SAC/{}/{}'.format(path,station,phase))
                    i +=1
                else:
                    pass
    else:
        print('The directory {}/Data/SAC_files/{} does not exists'.format(path,station))
        pass

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
#       Also lets load the gcarc from each stream, so we can test for whether SKKS should be measuable
        self.gcarc = (st[0].stats.sac.gcarc)
#       As this gcarc is calculated in split_read.py I know that it should be the same for all three traces
#       So for ease we will always read it from st[0]

    def model_traveltimes(self,phase):
        """
        Function to run TauP traveltime models for the SKS phase.
        Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object

        tr - trace object for which SKS arrival time will be predicted
        """
        model = ob.taup.tau.TauPyModel(model="iasp91")
        traveltime = model.get_travel_times(self.BHN[0].stats.sac.evdp,self.BHN[0].stats.sac.gcarc,[phase])[0].time

        return traveltime

    def check_phase_dist(self,phase):
        """
        Function to test if the given phase is actually measureable!
        """
        if phase == 'SKS':
            return True
        elif phase == 'SKKS':
            if self.gcarc >= 105.0:
                return True
            else:
                print('Event-Station distance less than 105 deg, too short for SKKS')
                return False
        else:
            print('Phase not SKS or SKKS')
            return False

    def process(self,station,phase,c1=0.01,c2=0.5):
        """
        Function to bandpass filter and trim the components
        t1 - [s] Lower bound of trim, time relative to event time
        t2 - [s] Upper bound of trim, time relative to event time
        c1 - [Hz] Lower corner frequency
        c2 - [Hz] Upper corner frequency
        By default traces will trim between 1300 - 1700s and will be filtered between 0.01Hz-0.5Hz
        """


#       De-mean and detrend each component
        # self.BHN.detrend(type='demean') #demeans the component
        # self.BHE.detrend(type='demean')
        # self.BHZ.detrend(type='demean')
#       Detrend
        # self.BHN.detrend(type='simple') #De-trends component
        # self.BHE.detrend(type='simple') #De-trends component
        # self.BHZ.detrend(type='simple') #De-trends component
#       Filter each component. Bandpass flag gives a bandpass-butterworth filter
        self.BHN.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHE.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHZ.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
#       Now trim each component to the input length
#       To ensure that we contain the phase information completlely lets model the arrival using TauP
        tt = self.model_traveltimes(phase)
#       Now set the trim
        t1 = (tt - 60) #I.e A minute before the arrival
        t2 = (tt + 120) #I.e Two minutes after the arrival
        self.BHN.trim(self.BHN[0].stats.starttime + t1,self.BHN[0].stats.starttime + t2)
        self.BHE.trim(self.BHE[0].stats.starttime + t1,self.BHE[0].stats.starttime + t2)
        self.BHZ.trim(self.BHZ[0].stats.starttime + t1,self.BHZ[0].stats.starttime + t2)
#       Add windowing ranges to sac headers user0,user1,user2,user3 [start1,start2,end1,end2]
#       Set the range of window starttime (user0/user1)
        user0 = tt - 15 #15 seconds before arrival
        user1 = tt # At predicted arrival
#       Set the raqnge of window endtime (user2/user3)
        user2 = tt + 15 # 15 seconds after, gives a min window size of 20 seconds
        user3 = tt + 30 # 30 seconds after, gives a max window size of 45 seconds
#
        self.BHN[0].stats.sac.user0,self.BHN[0].stats.sac.user1,self.BHN[0].stats.sac.user2,self.BHN[0].stats.sac.user3 = (user0,user1,user2,user3)
        self.BHE[0].stats.sac.user0,self.BHE[0].stats.sac.user1,self.BHE[0].stats.sac.user2,self.BHE[0].stats.sac.user3 = (user0,user1,user2,user3)
        self.BHZ[0].stats.sac.user0,self.BHZ[0].stats.sac.user1,self.BHZ[0].stats.sac.user2,self.BHZ[0].stats.sac.user3 = (user0,user1,user2,user3)

    def write_out(self,station,phase,i=0,path=None):
#       Now write out the three processed components
#       Naming depends on whether this is being executed as a test or within a loop
#       where a counter should be provided to prevent overwriting.

        if path is not None:
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

    def sheba(self,station,phase,i = 0,path=None):
        """
        The big one! This function uses the subprocess module to host sac and then runs sheba as a SAC macro
        """
        print('Passing {}_{}_{} into Sheba'.format(station,phase,i))

        p = sub.Popen(['sac'],
                     stdout = sub.PIPE,
                     stdin  = sub.PIPE,
                     stderr = sub.STDOUT,
                     cwd = path,
                     encoding='utf8')
#       Now that the child process SAC has been opened, lets specify what arguements we want to pass to it
#       echo on means commands should be echoed (so we can check whats going on if we print output)
#       SETMACRO makes sure SAC can find sheba (avoids pathing problems)
#       m sheba calls sheba as a SAC macro
        s = '''
        echo on\n
        SETMACRO /Users/ja17375/Ext_programs/macros
        m sheba file {}_{}_{} plot yes pick no nwind 10 10 batch yes
        '''.format(station,phase,i)
        try:
            out = p.communicate(s)
            # print(out[0])
        except CalledProcessError as err:
            print(err)


def tidy(station,phase):
    """
    Function to tidy up the working directory after a sheba run. Do things like concatenate final result files together, move postscripts to the postscript folder etc.
    Give the working directory a good clean basically
    """

    sub.call(['cat','', '{}_{}*.final_result'.format(station,phase)])
