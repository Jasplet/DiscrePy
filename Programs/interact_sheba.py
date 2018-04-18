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
import pandas as pd
import subprocess as sub
from subprocess import CalledProcessError
import os.path
import time
import shlex
from multiprocessing import Pool, current_process
import contextlib
from glob import glob
##############################
# Begin Program
##############################
def main(phase='SKS',batch=False,evt_sta_list=None):
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
        print(__name__)
        if __name__ == 'interact_sheba':
            ## Set of pool for mapping
            print('hi')
            statlist ='{}/Data/{}'.format(path,evt_sta_list)
            print('Processing Data from the Event-Station List {}'.format(statlist))
            stations = pd.read_csv(statlist,delim_whitespace=True).STAT.unique()

            out_pre = input('Enter SDB file name: ')
            outfile = '{}/Sheba/Results/{}_{}_sheba_results.sdb'.format(path,out_pre,phase)
            with contextlib.closing( Pool() ) as pool:
#           Iterate over stations in the station list.

                pool.map(run_sheba,stations)
#               pool.map(tidyup,stations) ??? Maybe this would work???
            tidy_path = 'Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Nulls_SKS'
            tidyup(tidy_path,phase,outfile)


    elif batch is False:
        station = input('Input Station Name > ')
        run_sheba(path,station,phase)


    end = time.time()
    runtime = end - start
    print('The runtime of main is {} seconds'.format(runtime))


def tidyup(path,phase,outfile):
    """
    Function to tidy up the working directory after a sheba run. Do things like concatenate final result files together, move postscripts to the postscript folder etc.
    Give the working directory a good clean basically
    """
 ## This needs to be reworked to be in python so that I can actually get this to work
    # sub.call(shlex.split('{}/Sheba/Programs/tidyup_by_stat.sh {} {} {} {}'.format(path,station,phase,path,outfile)))

    fnames = glob('{}/*/{}/*final_result'.format(path,phase))
    print(fnames)
    results = []
    for file in fnames:
        with open(file,'r') as input:
             head = input.readline()
             h = head.split()
             h.remove('%')
             i = h.index('STAT')
             h[2],h[3:i+1]= h[i],h[2:i]
             header = ' '.join(h)


             for line in input.readlines():
                r = line.split()
                r.remove('%')
                r[2],r[3:i+1] = r[i],r[2:i]
                result = ' '.join(r)
                results.append(result)

    results.insert(0,header)
    print('Writing Results to {}.sdb in /Users/ja17375/Shear_Wave_Splitting/Sheba/Results'.format(outfile))
    with open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}.sdb'.format(outfile),'w') as writer:
        for r in results:
            writer.write(str(r) + '\n')

def tidyup_final(path,phase,outfile):
    """
    Calls the bash script tidyup.sh to compile the by station results files into overall results.
    """
    ## This needs to be reworked in Python (get rid of call to bash function) so this will work in the new, parallel, framework!!!
    sub.call(shlex.split('{}/Sheba/Programs/tidyup.sh {} {}'.format(path,phase,outfile)))

def run_sheba(station,path='/Users/ja17375/Shear_Wave_Splitting',phase='SKS',outfile='null_results'):
    """
    Function that holds the guts of the workflow for preparing SAC files and running sheba
    """
    print('Worker {} running sheba for station {}'.format(current_process().pid,station))
    #Each station SHOULD have its own directory within Data/SAC_files
    #If the data has been downloaded. So lets look for directorys that exist


    dir_path = '{}/Data/SAC_files/{}'.format(path,station)
    if os.path.isdir(dir_path):
        #'Happy Days! The data directory exists!'
        i = 1 # Counter
        with open('{}/{}_downloaded_streams_Jacks_Nulls.txt'.format(dir_path,station),'r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
            for line in reader.readlines():
                line.strip('\n')
                label = line[55:-1].strip('/') # Extract the event label STAT_DATE_TIME so I can use it to label output stremas from sheba
                print(label)
                st_id = '{}BH?.sac'.format(str(line[0:-1]))
                st = ob.read(st_id)

                if len(st) is 3:
                    Event = Interface(st,station)
                    if Event.check_phase_dist(phase_to_check=phase) is True:
                        Event.process(station,phase)
                        outdir = '{}/Sheba/Runs/Jacks_Nulls_SKS/{}/{}'.format(path,station,phase)

                        try:

                            Event.write_out(phase,label,path=outdir)
                        except OSError:
                            print('Directory {} writing outputs do not all exist. Initialising'.format(outdir))
                            os.makedirs(outdir)
                            Event.write_out(phase,label,path=outdir)

                        Event.sheba(station,phase,label,i,path=outdir)

                        #tidyup_by_stat(path,station,phase,label,outfile)

                    else:
                        pass
#                   Counter, i , for number of events processed
                    i +=1
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
                print('Event-Station distance less than 105 deg, too short for SKKS')
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
        print('Passing {} ({}) into Sheba'.format(label,i))

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
