#! /usr/bin/env python
##############################
#   Program: sheba.py
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
#   This program is designed to be run from the command line (using python sheba.py)
#   If you want to run sheba for single stations or interactivly use parts of this program, then use the "module" version interact_sheba.py
#
#
##############################
#   Import Statements
##############################
#   Standard Packages - all freely available
import sys
import obspy as ob
import obspy.taup
import pandas as pd
import subprocess as sub
from subprocess import CalledProcessError
import os.path
import time
import shlex
from multiprocessing import Pool, current_process
import contextlib
from glob import glob
############################################################################################
# Define Functions and Classes
############################################################################################
def tidyup(path,phase,outfile):
    """
    Function to collect .final_result files output from Sheba into the Run directory. Results are written in the SBD format to Sheba/Results

    path - [str]: path to the Run directory (e.g. /users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Run_Name)

    phase - Phase that you want to collect results for

    outfile - [str]: outfile name (including full path)
    """
    fnames = glob('{}/*/{}/*final_result'.format(path,phase))
    results = []
    for i,file in enumerate(fnames):
        f_stat = fnames[i].rstrip('final_result') + 'stats'
        with open(file,'r') as input, open(f_stat,'r') as stats:
             head = input.readline()
             head_s = stats.readline()
             h = head.split()
             s = head_s.split()
             h.remove('%')
             s.remove('%')
             del s[0:2] ,s[-1]
             j = h.index('STAT')
             h[2],h[3:j+1]= h[j],h[2:j]
             header = ' '.join(h) + ' ' + ' '.join(s)
             print(header)

             for line in input.readlines():
                s = stats.readline().split() # Read the next line from the Stats file (should have the same number of lines as final_result)
                r = line.split()
                r.remove('%')
                s.remove('%')
                del s[0:2] ,s[-1]
                r[2],r[3:j+1] = r[j],r[2:j]
                result = ' '.join(r) + ' ' + ' '.join(s)
                results.append(result)

    results.insert(0,header)
    print('Writing Results to {}.sdb in /Users/ja17375/Shear_Wave_Splitting/Sheba/Results'.format(outfile))
    with open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}.sdb'.format(outfile),'w') as writer:
        for r in results:
            writer.write(str(r) + '\n')

def run_sheba(station,path='/Users/ja17375/Shear_Wave_Splitting',phases=['SKS','SKKS'],outfile='split_results'):
    """
    Function that holds the guts of the workflow for preparing SAC files and running sheba
    """
    #Each station SHOULD have its own directory within Data/SAC_files
    #If the data has been downloaded. So lets look for directorys that exist
    dir_path = '{}/Data/SAC_files/{}'.format(path,station)
    if os.path.isdir(dir_path):
        #'Happy Days! The data directory exists!'

        with open('{}/{}_downloaded_streams_Jacks_Split.txt'.format(dir_path,station),'r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
            for line in reader.readlines():
                for phase in phases:
                    #print(phase)
                    line.strip('\n')
                    label = line[55:-1].strip('/') # Extract the event label STAT_DATE_TIME so I can use it to label output stremas from sheba
                    st_id = '{}BH?.sac'.format(str(line[0:-1]))
                    st = ob.read(st_id)
                    if len(st) is 3:
                        Event = Interface(st,station)
                        if Event.check_phase_dist(phase_to_check=phase) is True:
                            Event.process(station,phase)
                            outdir = '{}/Sheba/Runs/Jacks_Split/{}/{}'.format(path,station,phase)
                            try:
                                Event.write_out(phase,label,path=outdir)
                            except OSError:
                                print('Directory {} writing outputs do not all exist. Initialising'.format(outdir))
                                os.makedirs(outdir)
                                Event.write_out(phase,label,path=outdir)

                            Event.sheba(station,phase,label,path=outdir)

                            #tidyup_by_stat(path,station,phase,label,outfile)

                        else:
                            pass

                    else:
                        print(' len(st) is not 3. Passing')
                        pass
    else:
        print('The directory {}/Data/SAC_files/{} does not exists'.format(path,station))
        pass

class Interface:
    """
    Class which will act as the interface to sheba.
    The "subprocess" sheba will be a bound method
    """
    def __init__(self,st,station):
        # self.date = date
        # self.time = time
        for i in [0,1,2]:
            st[i].stats.sac.kstnm = '{:>8}'.format(st[i].stats.sac.kstnm)
#§          Formats Station name in headers so that it is 8 characters long, with emtpy character fill with whitespaces
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

        self.station = station
#       As this gcarc is calculated in split_read.py I know that it should be the same for all three traces
#       So for ease we will always read it from st[0]

    def model_traveltimes(self,phase):
        """
        Function to run TauP traveltime models for the SKS phase.
        Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object

        tr - trace object for which SKS arrival time will be predicted
        """
        model = ob.taup.tau.TauPyModel(model="iasp91")
        # Add in a test for the case where depth has been gien in meters (as OBSPY IS DUMB AS HELL AND RETURNS DEPTH IN [m] FFS)
        if self.BHN[0].stats.sac.evdp >= 1000:
            #This might be a bit generous but my events shouldnt be shallower the ~ 1 km or deeper than 1000km anyway (and if this is the case than there is something SERIOUSLY wrong with our earth models)


            traveltime = model.get_travel_times((self.BHN[0].stats.sac.evdp/1000),self.BHN[0].stats.sac.gcarc,[phase])[0].time
        elif self.BHN[0].stats.sac.evdp == 0: # Theres an event where the event data couldnt be found so evdp was set to be 0
            # Having a depth of zero will give us problems so NOW change it to 10.0km exactly (these traveltimes could be very dodgy)
            err_out = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Events_with_evdp_of_0.txt','w+')
            err_out.write('Station: {}, has event starting at {} with an evdp of 0!\n'.format(self.station,self.BHN[0].stats.starttime))

            traveltime = model.get_travel_times(10,self.BHN[0].stats.sac.gcarc,[phase])[0].time

        else:
            traveltime = model.get_travel_times((self.BHN[0].stats.sac.evdp),self.BHN[0].stats.sac.gcarc,[phase])[0].time

        return traveltime

    def check_phase_dist(self,phase_to_check):
        """
        Function to test if the given phase is actually measureable!
        """
        if phase_to_check == 'SKS':
            return True
        elif phase_to_check == 'SKKS':
            if self.gcarc >= 105.0:
                return True
            else:
                #print('Event-Station distance less than 105 deg, too short for SKKS')
                return False
        else:
            print('Phase not SKS or SKKS')
            return False

    def process(self,station,phase,c1=0.01,c2=0.5):
        """
        Function to bandpass filter and trim the components
        Seismograms are trimmed so that they start 1 minute before the expected arrival and end 2 minutes after the arrival
        c1 - [Hz] Lower corner frequency
        c2 - [Hz] Upper corner frequency
        By default traces will be filtered between 0.01Hz-0.5Hz
        """


#       De-mean and detrend each component
#         self.BHN.detrend(type='demean') #demeans the component
#         self.BHE.detrend(type='demean')
#         self.BHZ.detrend(type='demean')
#         print(self.BHN)
# #       Detrend
#         self.BHN.detrend(type='simple') #De-trends component
#         self.BHE.detrend(type='simple') #De-trends component
#         self.BHZ.detrend(type='simple') #De-trends component
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


    def write_out(self,phase,label,path=None):
        """
        Function to write the component seismograms to SAC files within the sheba directory structure so Sheba can access them
        station [str] - station code
        phase [str] - phase code for which seismic phase splitting is being meausured
        i [int] - counter for nukber of events (with the desired phase) at the station
        path [str] - path that you want the seismogrmas saved to.
        """
#       Now write out the three processed components
#       Naming depends on whether this is being executed as a test or within a loop
#       where a counter should be provided to prevent overwriting.

        if path is not None:
            self.BHN.write('{}/{}{}.BHN'.format(path,label,phase),format='SAC',byteorder=1)
            self.BHE.write('{}/{}{}.BHE'.format(path,label,phase),format='SAC',byteorder=1)
            self.BHZ.write('{}/{}{}.BHZ'.format(path,label,phase),format='SAC',byteorder=1)
        else:
            self.BHN.write('{}.BHN'.format(label),format='SAC',byteorder=1)
            self.BHE.write('{}.BHE'.format(label),format='SAC',byteorder=1)
            self.BHZ.write('{}.BHZ'.format(label),format='SAC',byteorder=1)




    def plot_comp(self):
        """
        Quick Function to plot component together on one seismogram
        """
        st = self.BHN + self.BHE + self.BHZ
        st.plot(type='relative')

    def sheba(self,station,phase,label,i = 0,path=None):
        """
        The big one! This function uses the subprocess module to host sac and then runs sheba as a SAC macro
        """
        print('Worker {} Passing {} into Sheba for {}. Time is {}'.format(current_process().pid,label,phase,time.time()))

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
        m sheba file {}{} plot yes pick no nwind 10 10 batch yes
        '''.format(label,phase)
        try:
            out = p.communicate(s)
            # print(out[0])
        except CalledProcessError as err:
            print(err)

##############################
# Begin Program
##############################
if __name__ == '__main__':
    print(__name__)
#  Main Level, calls functions and set everything up
#  This function depends on the existence of a station list file specified by statlist and you have sac data alreayd downloaded
#  Path should point to the directory that you want the sheba processing directories to be stored under
#
#  Phases - [list]: Phases that you want to process for shear wave splitting using sheba. Traces should contain the expected arrivals. Interact sheb will loop over all phases provided and run sheba accordingly
#
#  Batch - [bool]: T- run in batch mode (multi station), using pythons multiprocessing library to improve performance. F - run for a single station
#
#  evt_sta_list - [str]: string pointing to the location of the event station list. This list should contain the a date/time and station that corresponds to the event time and station where arrivals are obseverd.
#                        This list is used to
#
# #First Indentify the possible station that we could have data for
#This way we know what directory paths to look in for the sac files

    #Loop over all stations in the list.
    #################### Time run #################
    start = time.time()
    #####################################################################################
    # Setup Paths, Constants and sort out input Arguements #############
    #####################################################################################
    # Set full path to station list
    evt_sta_list = sys.argv[1]

    statlist ='/Users/ja17375/Shear_Wave_Splitting/Data/{}'.format(evt_sta_list)
    # echo out where I expect the staiton list to be

    print('Processing Data from the Event-Station List {}'.format(statlist))
    stations = pd.read_csv(statlist,delim_whitespace=True).STAT.unique()
    # Get the user to input the outfile name they want (Phase label will be added later)
    out_pre = input('Enter SDB file name: ')
    ################### SET THE PHASES TO BE PROCESSED HERE #############################
    phases = ['SKS','SKKS']

    ######################################################################################
    ############### Run Sheba - using parallel processing with Pool ######################
    ######################################################################################
    with contextlib.closing( Pool(processes = 7) ) as pool:
    #           Iterate over stations in the station list.
        pool.map(run_sheba,stations)
    #               pool.map(tidyup,stations) ??? Maybe this would work???
    for phase in phases:
        """ Loop over phases process and tidyup results """
        tidy_path = 'Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split'
        outfile = '{}_{}_sheba_results.sdb'.format(phase,out_pre)
        tidyup(tidy_path,phase,outfile)

    ######################################################################################
    # End Timing of run
    ######################################################################################
    end = time.time()
    runtime = end - start
    print('The runtime of main is {} seconds'.format(runtime))

# END
