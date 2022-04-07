#------------------------------
# measure
#------------------------------
# This is a new and improved version of Split_Measure. Measuring capabilities are extend to Transverse Minimisation and Cross Correlation (along with Eigenvalue Minimisation).
# Each trace is windowed once, with this window being used for all 3 methods.
#
import obspy
import splitwavepy as sw
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import obspy.signal.filter as filt
import os.path
from obspy.clients.fdsn import Client
from obspy.taup import TauPyModel
import matplotlib.gridspec as gridspec
import Split_Measure as sm

def st_read(st_id):
    """
    Function to read sac files based on an input of expected file names (minus the compentend extension).
    It is expected that events will have a North, East and Vertical (BHN,BHE,BHZ) component. If there are not 3 files present then the event is treated as not being read properly.
    Events that are not read properly or that return an exception return a value of False to indicate that there is not sufficient data.
    """
    st = obspy.core.stream.read(st_id) # Reads all traces that match st_id
    if len(st) == 3: #Check is there are 3 traces in the stream (East,North and Vertical)
        return st
    else:
        return False

def read_pair(id):
    """
    Function to read a Pair/EigM file.
    Returns Pair object
    """
    data = sw.load(id)

    #Check if we are reading an EigenM object or a Pair
    if type(data) is 'splitwavepy.measure.eigenM.EigenM':
        quality = data.quality
        pair = data.data
        return pair, quality
    else:
        Exception('Unexpected Object type')

def TauP(source_depth_km, distance_degrees, phase_list=['SKS','SKKS']):
    """
    Function that uses the TauP toolkit to calculate traveltime for reaquired phases (default is that SKS and SKKS are calculated)
    """
    model = obspy.taup.tau.TauPyModel(model="iasp91") # Makes a TauPyModel class object which we use to model the traveltimes
    # Make list of traveltimes for input phases
    traveltimes = [model.get_travel_times(source_depth_km, distance_degrees,phase_list=phase_list)[0].time for phase in phase_list]

    return traveltimes

def process_st(st,tt,trim=120):
    """
    Function that filters, trims (so that there are an even number of points) and windows the traces. If pair does not exist it is created here
    Traces are filtered between 2 and 100 seconds (0.01Hz and 0.5Hz)
    """

    st.filter("bandpass",freqmin= 0.01, freqmax= 0.5,corners=2,zerophase=True)
    tt_UTC = st[0].stats.starttime + tt # phase predicted arrival as a UTCDateTimeobject so we can have it as the center of the trim.
    st.trim(tt_UTC-trim,tt_UTC + trim)

    pair = sw.Pair(st[1].data,st[0].data,delta= st[0].stats.delta)

    pair.plot(pick=True, marker = trim)
    w1,w2 = pair.wbeg(),pair.wend()
    wbeg_origin_time,wend_origin_time = (tt - trim + w1),(tt - trim + w2) # This should make my wbeg,wend relative to the origin time of the event)

    return pair, [wbeg_origin_time, wend_origin_time]

def splitting(pair,baz):
    """
    Function for carrying out the shear wave splitting analysis.
    All three methods from SplitwavePy (EigenM,TransM and CrossM) are performed here
    """
###################################################
#   Measure Splitiing using Eigenvalue Minimisation
    eigm = sw.EigenM(pair)
###################################################
#   Measure Splitting using Cross Correlation
    xcor = sw.CrossM(pair)
###################################################
#   Measure Splitting using transverse Minimisation
    trans = sw.TransM(pair,baz)

    return eigm, xcor, trans


def write(outfile,station,window,eigm=None,st=None,date=None,time=None,quality):

    if eigm is not None:
        #Splitting measure has been made or already exists and needs to be written out
        wl_results = sm.split_match(date,time,station)

        eig_results = [eigm.fast, eigm.dfast, eigm.lag, eigm.dlag] #Measurement that I want to output
        xcor_results = [xcor.fast,xcor.dfast,xcor.lag,xcor.dlag]
        trans_results = [trans.fast,trans.dfast,trans.lag,trans.dlag]

        attrib = ['stla','stlo','evla','evlo','evdp','gcarc','baz'] #SAC attribute values that I want to extract and print later
        stats = [st[0].stats.sac[i] for i in attrib] # Use list comprehension to extract sac attributes I want.
        org = [st[0].stats.station,date,time]

        outfile.write('{} {:07d} {:06d} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:06.2f} {:06.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:5.3f} {:4.2f} {} {} {} {} {} {} {}\n'.format(*org,*stats,*meas,str(eigm.quality[0])))
    elif eigm is None:
        meas = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN']
        stats = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN']
        org = [station,'NaN','NaN']
        quality = ['x']
        print('No stream for event')
        outfile.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22}\n'.format(*org,*stats,*meas,quality[0]))

def measure(line,phases=['SKS','SKKS']):

    st = st_read(st_id)

    if st!= False:
        origin_time = st[0].stats.starttime
#       Create time and date stamps
        date,time = int(str(origin_time.year)+str(origin_time.julday).zfill(3)),int(str(origin_time.hour).zfill(2)+str(origin_time.minute).zfill(2)+str(origin_time.second).zfill(2))

        eig_file ='{}{}_{:07d}_{:06d}.eigm'.format('/Users/ja17375/Python/SKS_Splitting/Eigm_Files/',station,date,time)
        if os.path.isfile(eig_file):
            pair = read_pair(eig_file)
        else:
            # Pair does not exist, so lets make one
            #for phase in phases: - Possible structuring for when I want to also consider SKKS
            traveltimes = TauP(source_depth_km=st[0].stats.sac.evdp,distance_degrees=st[0].stats.sac.evdp)

            pair, window = process_st(st,traveltimes[0])
        ####

        eigm, xcor, transm = splitting(pair,st[0].stats.sac.baz)

        write(outfile,station,date,time,st,window,eigm,xcor,transm,quality)
##########################
#   Past here all functions/classes are work in progress and should not be used!
def stack():
    """
    Function for stacking of error surfaces for different splitting measurements.
    """

class Quality():

    def __init__(self,*args,**kwargs):
