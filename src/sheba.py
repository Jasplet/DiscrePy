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
import os.path
import obspy as ob
import obspy.taup
import subprocess as sub
from subprocess import CalledProcessError
import time
import timeit
from multiprocessing import Pool, current_process
from functools import partial
import contextlib
from glob import glob
from Windower import WindowPicker
import argparse
############################################################################################
# Define Functions and Classes
############################################################################################
def tidyup(path,phase,outfile,outdir):
    """
    Function to collect .final_result files output from Sheba into the Run directory. Results are written in the SBD format to Sheba/Results

    path - [str]: path to the Run directory (e.g. /users/ja17375/DiscrePy/Sheba/Runs/Run_Name)
    phase - Phase that you want to collect results for
    outfile - [str]: outfile name (filename only! Path taken from Path variable)
    """
    if phase == 'SYNTH':
        fnames = glob('{}/*.final_result'.format(path))
        print('{}*.final_result'.format(path))
    else:
        fnames = glob('{}/*/{}/*final_result'.format(path,phase))
    results = []
    # print(fnames)
    for i,file in enumerate(fnames):
        print(file)
        f_stat = fnames[i].rstrip('final_result') + 'stats'

        with open(file,'r') as input, open(f_stat,'r') as stats:
             head = input.readline()
             head_s = stats.readline()
             h = head.split()
             s = head_s.split()
             # print(s)
             h.remove('%')
             s.remove('%')
             del s[0:2] ,s[-1]
             j = h.index('STAT')
             h[2],h[3:j+1]= h[j],h[2:j]

             header = ' '.join(h) + ' ' + ' '.join(s)

             for line in input.readlines():
                s = stats.readline().split() # Read the next line from the Stats file (should have the same number of lines as final_result)
                r = line.split()
                # if i == 0:
                    #'For the first file print out r
                    # print(r)
                r.remove('%')
                s.remove('%')
                del s[0:2] ,s[-1]
                r[2],r[3:j+1] = r[j],r[2:j]
                result = ' '.join(r) + ' ' + ' '.join(s)
                results.append(result)

    results.insert(0,header)
    # outdir = path.split('/')[-1]
    print('Writing Results to {} in /Users/ja17375/DiscrePy/Sheba/Results/{}'.format(outfile,outdir))
    with open('/Users/ja17375/DiscrePy/Sheba/Results/{}/{}'.format(outdir,outfile),'w') as writer:
        for r in results:
            writer.write(str(r) + '\n')

def run_synth_i(runpath,filepath):
    '''
    Self Contained Function to run sheba (DEISGNED FOR IPYTHON USE) for a set of synthetics generated by sacsplitwave
    '''
    r_dir = '/Users/ja17375/DiscrePy/Sheba/Runs/SYNTH/LowSNR'
    with open('/Users/ja17375/DiscrePy/Data/Synthetics.events','r') as reader:
        for i, line in enumerate(reader.readlines()):
            file = line.strip('\n')

            label = '{}_'.format(file.split('/')[-1])
            # print(label)
            st_id = '{}.BH?'.format(file)
            st = obspy.read(st_id)
            stat= st[0].stats.station
            #print(st)
            Event = Interface(st)
            Event.process(synth=True)
            # print(r_dir)
            print(label)
            Event.write_out('SYNTH',label,path=r_dir)
            Event.sheba(stat,'SYNTH',label,path=r_dir)
     #End of loop
    tidyup(r_dir,'SYNTH','Synthetics_results.sdb')

def run_synth(runpath,filepath):
    '''
    Runs sheba for a list of file mapped in (for parallelised version).
    '''
    # print('Running')
    label = '{}_'.format(filepath.split('/')[-1])
    phase= 'SYNTH' #Empty String for phase are these are synthetics
    # print(label)
    st_id = '{}.BH?'.format(filepath)
    print(st_id)
    st = obspy.read(st_id)
    #print(st)
    Event = Interface(st)
    # print('Synthetics Used, Windows *should* be predefined')
    Event.process(synth=True)
    # print(r_dir)
    print(label)
    Event.write_out('SYNTH',label,path=runpath)
    # print('{}, {}, {}, {}, False'.format(stat,phase,label,runpath))
    Event.sheba(phase,label,path=runpath,nwind=True)



def run_sheba(runpath,filepath,phases=['ScS']):
    """
    Function that holds the guts of the workflow for preparing SAC files and running sheba
    """
    #Each station SHOULD have its own directory within Data/SAC_files
    #If the data has been downloaded. So lets look for directorys that exist

    dir_path = '/'.join(filepath.split('/')[0:-1])
    # print('RP {} \n FP: {} '.format(runpath,filepath))
    if os.path.isdir(dir_path):
        #'Happy Days! The data directory exists!'
        for phase in phases:
            #print(phase)
            label = '{}'.format(filepath.split('/')[-1]) # Extract the event label STAT_DATE_TIME so I can use it to label output stremas from sheba
            # print('Label is {}'.format(label))
            st_id = '{}BH?.sac'.format(filepath)
            # st_id = '{}.BX?'.format(filepath) # .BX? for SPECFeM SYNTHETICS
            st = ob.read(st_id)
            station = st[0].stats.station
            f_check = '{}/{}/{}/{}{}_sheba.final_result'.format(runpath,station,phase,label,phase)
            # print('Fcheck is {}'.format(f_check))
            # print('{}/{}/{}/{}{}.B?E'.format(runpath,station,phase,label,phase))
            if os.path.isfile(f_check) == True: # Check if event has already been processed
                print('File has already been processed: {} '.format(f_check))
            else:
                #print('File to process: {} '.format(f_check))
                if len(st) == 3:
                    Event = Interface(st)
                    print('Check epicentral distances')
                    if Event.check_phase_dist(phase_to_check=phase) is True:
                #       To ensure that we contain the phase information completlely lets model the arrival using TauP
                        print('Modelling Traveltimes')
                        Event.model_traveltimes(phase)
                        if Event.tt is None:
                            print('NO Traveltimes, bad data file {}'.format(st_id))
                            pass
                        else:
                            # print('Process waveform, stat {}, phase {}, label {}'.format(station,phase,label))
                            Event.process(synth=False,window=False) # If windowing is set to true then SHEBA must be called in serial mode
                            if Event.bad is True:
                                # A bad waveform that we dont want anyhting to do with, so skip it
                                pass

                            else:
                                # print('Function', Event.process(phase,synth=False))
                                outdir = '{}/{}/{}'.format(runpath,station,phase)
                                try:
                                    Event.write_out(phase,label,path=outdir)
                                except OSError:
                                    print('Directory {} writing outputs do not all exist. Initialising'.format(outdir))
                                    os.makedirs(outdir)
                                    # print('Label is {}. Path is {}'.format(label,path))
                                    Event.write_out(phase,label,path=outdir)
                                # print('RUn sheba, stat {}, phase {}, label {}, out {}'.format(station,phase,label,outdir))
                                Event.sheba(phase,label,path=outdir,nwind=True)
                                #tidyup_by_stat(path,station,phase,label,outfile)
                    else:
                        # print('Fail, Distance ')
                        pass
                else:
                    print(' len(st) is not 3. Passing')
                    pass
    else:
        print('The directory {} does not exists'.format(dir_path))
        pass

class Interface:
    """
    Class which will act as the interface to sheba.
    The "subprocess" sheba will be a bound method
    """
    def __init__(self,st):
        # self.date = date
        # self.time = time
        self.ch = 'H'
        for i in [0,1,2]:
            st[i].stats.sac.kstnm = '{:>8}'.format(st[i].stats.sac.kstnm)
#          Formats Station name in headers so that it is 8 characters long, with emtpy character fill with whitespaces
        try:
            self.BHE = st.select(channel='BHE')
            self.BHE[0].stats.sac.cmpinc = 90
            self.BHE[0].stats.sac.cmpaz = 90
            self.BHN = st.select(channel='BHN')
            self.BHN[0].stats.sac.cmpinc = 90
            self.BHN[0].stats.sac.cmpaz = 0
            self.BHZ = st.select(channel='BHZ')
            self.BHZ[0].stats.sac.cmpinc = 0
            self.BHZ[0].stats.sac.cmpaz = 0

        except IndexError:
            try:
                print('BH? channels not found, trying BX?')
                self.BHE = st.select(channel='BXE')
                self.BHE[0].stats.sac.cmpinc = 90
                self.BHE[0].stats.sac.cmpaz = 90
                self.BHN = st.select(channel='BXN')
                self.BHN[0].stats.sac.cmpinc = 90
                self.BHN[0].stats.sac.cmpaz = 0
                self.BHZ = st.select(channel='BXZ')
                self.BHZ[0].stats.sac.cmpinc = 0
                self.BHZ[0].stats.sac.cmpaz = 0
                self.ch = 'X'
            except IndexError:
                    print('Channel selection failed, indexing instead')
                    self.BHE = st[0]
                    self.BHE.stats.sac.cmpinc = 90
                    self.BHE.stats.sac.cmpaz = 90
                    self.BHN = st[1]
                    self.BHN.stats.sac.cmpinc = 90
                    self.BHN.stats.sac.cmpaz = 0
                    self.BHZ = st[2]
                    self.BHZ.stats.sac.cmpinc = 0
                    self.BHZ.stats.sac.cmpaz = 0
        # print(self.BHN)
#       Also lets load the gcarc from each stream, so we can test for whether SKKS should be measuable
        print('Checking for lost keys')
        if all (k in self.BHE[0].stats.sac for k in ('user0','user1','user2','user3')):
            print('No lost keys!')
        else:
            print('We\'ve lost at least one key, added to dict')
            keychain = {'user0':0,'user1':1,'user2':2,'user3':3}
            self.BHE[0].stats.sac.update(keychain)
            self.BHN[0].stats.sac.update(keychain)
            self.BHZ[0].stats.sac.update(keychain)
        self.gcarc = (st[0].stats.sac.gcarc)
        self.station = st[0].stats.station
        self.delta = st[0].stats.delta
        self.bad = False
        # self.path = os.getcwd()
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
            print(traveltime)
        elif self.BHN[0].stats.sac.evdp == 0: # Theres an event where the event data couldnt be found so evdp was set to be 0
            # Having a depth of zero will give us problems so NOW change it to 10.0km exactly (these traveltimes could be very dodgy)
            err_out = open('/Users/ja17375/DiscrePy/Sheba/Events_with_evdp_of_0.txt','w+')
            err_out.write('Station: {}, has event starting at {} with an evdp of 0!\n'.format(self.station,self.BHN[0].stats.starttime))
            traveltime = model.get_travel_times(10,self.BHN[0].stats.sac.gcarc,[phase])[0].time
        else:
            tt = model.get_travel_times((self.BHN[0].stats.sac.evdp),self.BHN[0].stats.sac.gcarc,[phase])
            # print(self.BHN[0].stats.sac)
            # print(tt)
            try:
                traveltime = tt[0].time
            except IndexError:
                print('index Error')
                traveltime =None

        evt_time = obspy.UTCDateTime(year = self.BHN[0].stats.sac.nzyear, julday = self.BHN[0].stats.sac.nzjday,hour=self.BHN[0].stats.sac.nzhour,minute=self.BHN[0].stats.sac.nzmin,second=self.BHN[0].stats.sac.nzsec,microsecond=self.BHN[0].stats.sac.nzmsec)
        self.tt = evt_time + traveltime
        self.tt_rel = traveltime
        return traveltime

    def check_phase_dist(self,phase_to_check):
        """
        Function to test if the given phase is actually measureable!
        """
        if phase_to_check == 'SKS':
            if self.gcarc < 145.0:
                return True
            else:
                print('Event-Station distance further than 145 deg, too far for SKS')
                return False

        elif phase_to_check == 'SKKS':
            if self.gcarc >= 105.0:
                return True
            else:
                print('Event-Station distance less than 105 deg, too short for SKKS')
                return False
        elif phase_to_check == 'ScS':
            if self.gcarc <= 95.0:
                return True
            else:
                print('Event-Station distance greater than 105 deg, too far for ScS')
        else:
            print('Phase {} not SKS or SKKS'.format(phase_to_check))
            return False

    def process(self,synth=False,c1=0.01,c2=0.5,window=False):
        """
        Function to bandpass filter and trim the components
        Seismograms are trimmed so that they start 1 minute before the expected arrival and end 2 minutes after the arrival
        c1 - [Hz] Lower corner frequency
        c2 - [Hz] Upper corner frequency (N.B I have used c2 = 00.5 intially and am now trying c2 = 0.3 to see if that improves SNR without cutting oiut too mcuh singal)
        By default traces will be filtered between 0.01Hz-0.5Hz
        """
        # print(synth)
      # De-mean and detrend each component
        self.BHN.detrend(type='demean') #demeans the component
        self.BHE.detrend(type='demean')
        self.BHZ.detrend(type='demean')
        # print(self.BHN)
#       Detrend
        self.BHN.detrend(type='simple') #De-trends component
        self.BHE.detrend(type='simple') #De-trends component
        self.BHZ.detrend(type='simple') #De-trends component
#       Filter each component. Bandpass flag gives a bandpass-butterworth filter
        self.BHN.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHE.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
        self.BHZ.filter("bandpass",freqmin= c1, freqmax= c2,corners=2,zerophase=True)
#       Now trim each component to the input length

        if synth == False: # We only need to trim and set the window length for real data, not synthetics made with sacsplitwave
#       Now set the trim

            print('Trim Traces')
            # print('tt = {}, Start = {}, End = {}'.format(self.tt,self.BHN[0].stats.starttime,self.BHN[0].stats.endtime ))
            t1 = (self.tt - 60) #I.e A minute before the arrival
            t2 = (self.tt+ 120) #I.e Two minutes after the arrival
            # print('t1 = {}, t2 = {}'.format(t1,t2))
            # self.BHN.trim(self.BHN[0].stats.starttime + t1,self.BHN[0].stats.starttime + t2)
            # self.BHE.trim(self.BHE[0].stats.starttime + t1,self.BHE[0].stats.starttime + t2)
            # self.BHZ.trim(self.BHZ[0].stats.starttime + t1,self.BHZ[0].stats.starttime + t2)
            self.BHN.trim(t1,t2)
            self.BHE.trim(t1,t2)
            self.BHZ.trim(t1,t2)
            # print('BX? records are assumed to already by short, so no trimming')
    #       Add windowing ranges to sac headers user0,user1,user2,user3 [start1,start2,end1,end2]
            # print(self.BHN[0].stats.sac.user0)
            if self.BHN[0].stats.sac.user0 == 0:
                print("Setting Window start/end ranges")
                # Set the range of window starttime (user0/user1)
                user0 = self.tt_rel - 15 #15 seconds before arrival
                user1 = self.tt_rel # t predicted arrival
        #       Set the raqnge of window endtime (user2/user3)
                user2 = self.tt_rel + 15 # 15 seconds after, gives a min window size of 20 seconds
                user3 = self.tt_rel + 30 # 30 seconds after, gives a max window size of 45 seconds
                keychain = {'user0':user0,'user1':user1,'user2':user2,'user3':user3}
                self.BHE[0].stats.sac.update(keychain)
                self.BHN[0].stats.sac.update(keychain)
                self.BHZ[0].stats.sac.update(keychain)
            else:
                print("Windows already set, user0-3 already set")
                print("User0 = ", self.BHN[0].stats.sac.user0)
                print("User1 = ",self.BHN[0].stats.sac.user1)
                print("User2 = ",self.BHN[0].stats.sac.user2)
                print("User3 = ",self.BHN[0].stats.sac.user3)

            if window == True:
                # Windowing code
                # Combine BHN and BHE to make a stream
                st = self.BHE + self.BHN
                # print(user0,user1,user2,user3)
                Windower = WindowPicker(st,user0,user1,user2,user3,self.tt_rel)
                # Windower.pick()
                if Windower.wbeg1 is None:
                    print("Skipping")
                    self.bad = True # Switch to tell sheba.py if we actually want to meausre this event
                else:
                    print("Windower Closed, adjusting window ranges")
                    (user0,user1,user2,user3) = Windower.wbeg1, Windower.wbeg2, Windower.wend1, Windower.wend2
                    self.bad = False
                # Set window ranges in SAC headers
                self.BHN[0].stats.sac.user0,self.BHN[0].stats.sac.user1,self.BHN[0].stats.sac.user2,self.BHN[0].stats.sac.user3 = (user0,user1,user2,user3)
                self.BHE[0].stats.sac.user0,self.BHE[0].stats.sac.user1,self.BHE[0].stats.sac.user2,self.BHE[0].stats.sac.user3 = (user0,user1,user2,user3)
                self.BHZ[0].stats.sac.user0,self.BHZ[0].stats.sac.user1,self.BHZ[0].stats.sac.user2,self.BHZ[0].stats.sac.user3 = (user0,user1,user2,user3)
        else:
            # print('Synthetics Used, Windows *should* be predefined')
            # print(self.BHN.stats.sac.user0,self.BHN.stats.sac.user1,self.BHN.stats.sac.user2,self.BHN.stats.sac.user3)
            pass

    def write_out(self,phase,label,path=None,synth=False):
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
            self.BHN.write('{}/{}{}.BHN'.format(path,label,phase,self.ch),format='SAC',byteorder=1)
            self.BHE.write('{}/{}{}.BHE'.format(path,label,phase,self.ch),format='SAC',byteorder=1)
            self.BHZ.write('{}/{}{}.BHZ'.format(path,label,phase,self.ch),format='SAC',byteorder=1)
        elif synth == True:
            self.BHN.write('{}/{}SYNTH.BHN'.format(path,label,self.ch),format='SAC',byteorder=1)
            self.BHE.write('{}/{}SYNTH.BHE'.format(path,label,self.ch),format='SAC',byteorder=1)
            self.BHZ.write('{}/{}SYNTH.BHZ'.format(path,label,self.ch),format='SAC',byteorder=1)
        else:
            self.BHN.write('{}.BHN'.format(label,self.ch),format='SAC',byteorder=1)
            self.BHE.write('{}.BHE'.format(label,self.ch),format='SAC',byteorder=1)
            self.BHZ.write('{}.BHZ'.format(label,self.ch),format='SAC',byteorder=1)

    def plot_comp(self):
        """
        Quick Function to plot component together on one seismogram
        """
        st = self.BHN + self.BHE + self.BHZ
        st.plot(type='relative')

    # def gen_infile(self,path,label,phase,nwind=10,tlag_max=4.0):
    #
    #     os.chdir(path) # Make sure we are in the right directory
    #
    #     with open('sheba.in','w') as writer:
    #         writer.write('SHEBA.IN \n')
    #         writer.write('{}{} \n'.format(label,phase)) # write file prefix
    #         writer.write('{} \n'.format(self.BHE[0].stats.channel)) # Write channels for each component (E, N, Z order)
    #         writer.write('{} \n'.format(self.BHN[0].stats.channel))
    #         writer.write('{} \n'.format(self.BHZ[0].stats.channel))
    #         writer.write('1 \n') # Specifies Eigenvalue minimisation, replace with spol if transverse minimisation is desired (not supported here)
    #         writer.write('{:i} {:i} \n'.format(nwind,nwind))
    #         writer.write('{} \n'.format(tlag_max)) # sets max tlag in gridsearch
    #         writer.write('0')

    def sheba(self,phase,label,path=None,nwind=True,):
        """
        The big one! This function uses the subprocess module to host sac and then runs sheba as a SAC macro
        """
        print('Worker {} Passing {} into Sheba for {}. Time is {}'.format(current_process().pid,label,phase,time.ctime()))
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
        if nwind == True:

            s = '''
            echo on\n
            SETMACRO /Users/ja17375/Ext_programs/macros
            m sheba file {}{} plot yes pick no nwind 10 10 batch yes
            '''.format(label,phase)
        elif nwind == False:
            s = '''
            echo on\n
            SETMACRO /Users/ja17375/Ext_programs/macros
            m sheba file {}{} plot yes pick yes batch yes
            '''.format(label,phase)
        try:
            # print(s)
            out = p.communicate(s)
            # print(out[0])
        except CalledProcessError as err:
            print(err)
##############################
# Begin Program
##############################
if __name__ == '__main__':
    # print(__name__)
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
    start = timeit.default_timer()
    #####################################################################################
    # Setup Paths, Constants and sort out input Arguements #############
    #####################################################################################
    # Set full path to station list
    parser = argparse.ArgumentParser()
    parser.add_argument("-e","--events",action="store",required=True,help="File name for the .events file containing the filestems you want to processes")
    parser.add_argument("-r","--rundir",action="store",required=True,help="The run directory that you wish to work in")
    parser.add_argument("-p","--parallel",action="store_true",default=False,help="Swtich for whether you want sheba to run in parallel")
    parser.add_argument("-nc","--ncores",action="store",type=int,help="Number of Cores requested for a parallel job")
    parser.add_argument("-s","--syn",action="store_true",default=False,help="Indicates that the data to be processed is synthetics, which can required a slightly different treatment")
    args = parser.parse_args()

    print(args)
    # evt_list = sys.argv[1] # A .events file containing a list of filestems to al the data we want to measure
    # rundir=sys.argv[2] # The run directory that you want to house the output files
    # run_mode = sys.argv[3] # par is wanting to run in parallel, ser if running serially
    # mode=sys.argv[4] # data if using real data, syn is using synthetics
    # print(mode,type(mode))
    file_list ='/Users/ja17375/DiscrePy/Data/{}'.format(args.events)
    # echo out where I expect the staiton list to be
    print('Processing Data from the Downloaded Event List {}'.format(file_list))
    files = []
    with open(file_list,'r') as reader:
        for i, line in enumerate(reader.readlines()):
            f = line.strip('\n')
            files.append(f)
            # print(f)
        print('There are {} files to process'.format(i+1))
    # Get the user to input the outfile name they want (Phase label will be added later)
    out_pre = input('Enter SDB file name: ')
    ################### SET THE PHASES TO BE PROCESSED HERE #############################
    # phases = ['SKS','SKKS']
    phases = ['ScS']
    ######################################################################################
    ############### Run Sheba - using parallel processing with Pool ######################
    ######################################################################################
    runpath ='/Users/ja17375/DiscrePy/Sheba/Runs/{}'.format(args.rundir)
    print('Runpath is :',runpath)
    if args.syn is False:
        # runpath ='/Users/ja17375/DiscrePy/Sheba/Runs/{}'.format(rundir)
        runner = partial(run_sheba,runpath,phases=['ScS'])
        if args.parallel is True:
            print("Parallel job. {} cores requested".format(args.ncores))
            with contextlib.closing( Pool(processes = args.ncores) ) as pool:
            #           Iterate over stations in the station list.
                pool.map(runner,files)
            #               pool.map(tidyup,stations) ??? Maybe this would work???
    #       Tidy up results
            print('Sheba run complete, time to tidy up')
        else:
            print("Running Serial Job")
            #Run in serial mode (booooo)
            for file in files:
                # print(file)
                run_sheba(runpath=runpath,filepath=file)

        for phase in phases:
            """ Loop over phases process and tidyup results """
            tidy_path = '/Users/ja17375/DiscrePy/Sheba/Runs/{}'.format(args.rundir)
            outfile = '{}_{}_sheba_results.sdb'.format(out_pre,phase)
            # outdir = tidy_path.split('/')[-1]
            tidyup(tidy_path,phase,outfile,args.rundir)


    elif args.syn is True:
        phase='SYNTH'
        runner = partial(run_synth,runpath)
        # As all the synthetics files are in the same directory (They share a "station") we cannot run in parralel with all the synthetics being stored in the same directory.
        # This occurs as my sythetics are all stored under one station directory, so we cannot keep each worker looking at a different station
        # As a result, the workers trip over each other an corrupt the majority of infiles.
        print('Synthetics Run')
        print('Runpath is :',runpath)
        # with contextlib.closing( Pool(processes = 8) ) as pool:
        #           Iterate over stations in the station list.
            # pool.map(runner,files)
        for file in files:
            runner(file)
        #Tidyup results
        print('Sheba run complete, time to tidy up')
        """ Loop over phases process and tidyup results """
        tidy_path = '/Users/ja17375/DiscrePy/Sheba/Runs/{}'.format(args.rundir)
        outfile = '{}_SYNTH_sheba_results.sdb'.format(out_pre)
        # outdir = rundir.split('/')[0]
        tidyup(tidy_path,phase,outfile,args.rundir)

    # print('Sheba run complete, time to tidy up')
    ######################################################################################
    ################ Run Sheba - Serial Model ########
    ######################################################################################
    # for i,file in enumerate(files):
    #     run_sheba(file,'/Users/ja17375/DiscrePy/Sheba/Runs/Test',phases= ['SKKS'])
    #     print('Iteration: {}. File: {}'.format(i,file))

    ######################################################################################
    # End Timing of run
    ######################################################################################
    end = timeit.default_timer()
    runtime = end - start
    print('The runtime of main is {} minutes'.format(runtime/60))
# END