#! /usr/bin/env python
## This is a pseudo package (a bundle of functions). Eventually this should turn into something more structured (eventually) such as a class?

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
import matplotlib.gridspec as gridspec
import os.path

###################################################

def splitting(station,switch,files,phase):
    """
    Measures SKS splitting for all streams listed in a ttext file at the provided path. These streams must be saved as SAC files.abs
    This function is the primary part of this module/package/script/thing, all the pther functions support this one.

    inpath - string contaiing the path of the textfile which contains the sac file to be read

    switch - optional kwarg to specify if you want to manually window the data or use a set of windows (Walpoles windows are availbale for use by default)
    """
    outfile = output_init(station,switch,phase)


    with open(files,'r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
        for line in reader.readlines():
            line.strip('\n')
            st_id = '{}*'.format(str(line[0:-1]))

            st = read_sac(st_id)
            # Intialise some global variables which I need to pass things between fucntions (this is probably not be best way to do things but it works!)
            global pair_glob
            global quality
            if st != False: # i.e. if the stream is sufficiently populated and has been read.

                SKS_UTC, t0, SKS = model_traveltimes(st[0],phase) # Returns SKS arrival as a UTCDateTime object, event origin time and SKS arrival relative to the origin time
                # print(t0)
                # print('SKS_UTC ={}'.format(SKS_UTC))
                quality = [] # variable to hold Callback key entries for estimated quality of splitting measurements
                date,time = int(str(t0.year)+str(t0.julday).zfill(3)),int(str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+str(t0.second).zfill(2)) #creates time and date stamps
                if switch is 'on':
                    eig_file ='{}/{}/{}/{}_{:07d}_{:06d}.eigm'.format('/Users/ja17375/Python/Shear_Wave_Splitting/Eigm_Files/',phase,station,station,date,time)
                elif switch is 'off':
                    eig_file = '{}/{}/{}/{}_{}_{:07d}_{:06d}'.format('/Users/ja17375/Python/Shear_Wave_Splitting/Eigm_Files/',phase,station,station,'JW_Windows',date,time)

                if os.path.isfile(eig_file):
                        print('Splitting already measured for this event, skipping')
                        split = sw.load(eig_file) #loads eigm file
                        # if split.quality == None:
                        # split.quality = 't'

                        write_splitting(outfile,station,eigm=split,st=st,date=date,time=time)

                else:
                    pair = st_prep(st = st, f_min = 0.01,f_max = 0.5)
                    # print(st[0].stats.starttime)
                    pair_glob = pair

                    if switch == 'on': # If manual windowing is on
                        split, wbeg, wend,fig = split_measure(pair,SKS)
                        split.quality = quality
                        plt.savefig('{}{}_{}_{:07d}_{:06d}'.format('/Users/ja17375/Python/Shear_Wave_Splitting/Figures/Eigm_Surface/',station,phase,date,time))
                        plt.close()

                    elif switch == 'off': #Manual windowing is off. For now this will just mean Jacks windows will be used. Eventually add automation or support for entering windows.

                        (wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend) = split_match(date,time,station)
                        pair.set_window(wl_wbeg,wl_wend) # Sets window to that of Jacks
                        wbeg,wend = pair.wbeg(), pair.wend()
                        split = sw.EigenM(pair,lags=(4,) ) # Measures splitting

                        fig = plt.figure(figsize=(12,6))
                        eigen_plot(split,fig)
                        plt.savefig('{}/{}_{}_{}_{:07d}_{:06d}'.format('/Users/ja17375/Python/Shear_Wave_Splitting/Figures/Eigm_Surface',station,phase,'JW_Windows',date,time))
                        plt.close()
                        split.quality = 'w'


                    write_splitting(outfile,station,eigm=split,st=st,date=date,time=time) #Call write_splitting where there are measuremtns to output
                    # save_sac(st,quality[0],date,time,wbeg,wend,switch)
                    split.save(eig_file) # Saves splitting measurements
            else:
                write_splitting(outfile,station)

    plt.close('all')
    outfile.close()

def write_splitting(outfile,station,eigm=None,st=None,date=None,time=None):

    if phase == 'SKS':

        if eigm is not None:
            #Splitting measure has been made or already exists and needs to be written out
            (wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend) = split_match(date,time,station)

            meas = [eigm.data.wbeg(), eigm.data.wend(), eigm.fast, eigm.dfast, eigm.lag, eigm.dlag,wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend ] #Measurement that I want to output
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

    elif phase == 'SKKS':

        if eigm is not None:
            meas = [eigm.data.wbeg(), eigm.data.wend(), eigm.fast, eigm.dfast, eigm.lag, eigm.dlag] #Measurement that I want to output
            attrib = ['stla','stlo','evla','evlo','evdp','gcarc','baz'] #SAC attribute values that I want to extract and print later
            stats = [st[0].stats.sac[i] for i in attrib] # Use list comprehension to extract sac attributes I want.
            org = [st[0].stats.station,date,time]
            outfile.write('{} {:07d} {:06d} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:06.2f} {:06.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:5.3f} {:4.2f} {}\n'.format(*org,*stats,*meas,str(eigm.quality[0])))
        elif eigm is None:
            meas = ['NaN','NaN','NaN','NaN','NaN','NaN']
            stats = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN']
            org = [station,'NaN','NaN']
            quality = ['x']
            print('No stream for event')
            outfile.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}\n'.format(*org,*stats,*meas,quality[0]))

def output_init(station,switch,phase):
    """
    Initialises output variables and output textfile

    station - string containing the station code
    """
    if switch is 'on':
        default_out = '/Users/ja17375/Python/Shear_Wave_Splitting/Measurements/{}_{}_Splitting.txt'.format(station,phase) #Default output filename
    elif switch is 'off':
        default_out = '/Users/ja17375/Python/Shear_Wave_Splitting/Measurements/{}_{}_Splitting_JW_Windows.txt'.format(station,phase)

    if os.path.isfile(default_out):
        #Default file exists! Request user permission to overwrite
        ovr = user_in('c1',default_out)
        if ovr == 'y':
            outfile = open(default_out,'w+')

        elif ovr == 'n':
            new_out = user_in('c2',station)
            outfile = open('/Users/ja17375/Python/Shear_Wave_Splitting/Measurements/{}.txt'.format(new_out),'w+')
    elif os.path.isfile(default_out) == False:
        print('{} does not exist and will be created'.format(default_out))
        outfile = open(default_out,'w+')

    if phase == 'SKS':
        outfile.write('STAT DATE TIME STLA STLO EVLA EVLO EVDP GCARC BAZ WBEG WEND FAST DFAST TLAG DTLAG WL_FAST WL_DFAST WL_TLAG WL_DTLAG WL_WBEG WL_WEND QUAL\n')
    elif phase == 'SKKS':
        outfile.write('STAT DATE TIME STLA STLO EVLA EVLO EVDP GCARC BAZ WBEG WEND FAST DFAST TLAG DTLAG QUAL\n')
    return outfile

def user_in(case,file1,station=None):
    if case == 'c1':
        a = input('The file {} already exists, do you want to overwrite? (y/n)\n'.format(file1))
        return a
    elif case == 'c2':
        b = input('Please input the desired file name (extension not required): \n')
        c = input('The filename will be {}.txt then? (y/n)'.format(b))
        if c == 'y':
            return b
        elif c == 'x':
            user_in(c2,station)

def read_sac(st_id):
    """
    Function to read sac files based on an input of expected file names (minus the compentend extension).
    It is expected that events will have a North, East and Vertical (BHN,BHE,BHZ) component. If there are not 3 files present then the event is treated as not being read properly.
    Events that are not read properly or that return an exception return a value of False to indicate that there is not sufficient data.
    """
    # try:
    # print(st_id)
    st = obspy.core.stream.read(st_id) # Reads all traces that match st_id
    if len(st) == 3: #Check is there are 3 traces in the stream (East,North and Vertical)
        return st
    else:
        return False
    # except Exception:
    #     print("Exception Encountered")
    #     return False


def save_sac(st,qual,date,time,wbeg,wend,switch):
    """
    Function to trim sac traces once they have been windowed and save these windowed traces for re-use. This is to add easy repeatability
    Arguements:
    st - the stream to process and save
    date - datestamp in the format yyyjjj (julian day)
    time - timestamp in format hhmmss
    wbeg - start of the window
    wend - end of the window
    qual - quaulitativr estimate of seismogram quilty (how clear the SKS arrival is)
    """
    if switch == 'on':
        for tr in st:
            ch = tr.stats.channel
            stat = tr.stats.station
            t0 = tr.stats.starttime
            path = '/Users/ja17375/Python/Shear_Wave_Splitting/Data/Proccessed_Streams'
            tr2 = tr.trim(t0 + wbeg, t0 + wend)
            tr2.write('{}/{}_{}_{:07d}_{:06d}_{}.sac'.format(path,stat,qual,date,time,ch), format='SAC')
    elif switch == 'off':
        pass

def model_traveltimes(tr,phase):
    """
    Function to run TauP traveltime models for the SKS phase.
    Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object

    tr - trace object for which SKS arrival time will be predicted
    """
    model = obspy.taup.tau.TauPyModel(model="iasp91")
    travelt = model.get_travel_times(tr.stats.sac.evdp,tr.stats.sac.gcarc,[phase])[0].time
    t0 = tr.stats.starttime # Start time of stream for event. This should be the event origin time.
    travelt_UTC = obspy.core.utcdatetime.UTCDateTime(t0 + travelt)# SKS arrival time relative to trace start as a UTCDateTime object

    return travelt_UTC, t0, travelt


def st_prep(st,f_min,f_max):
    """
    Prepares Stream for spltting analysis (by bandpass filtering and trimming) and then creates the Pair object
    """

    if st[0].data.size != st[1].data.size: #Check that the streams are the same length
        len_diff = abs(st[0].data.size - st[1].data.size) # length difference in number of point
        if  st[0].data.size  >  st[1].data.size : # st[0] longer so trim from its end
            st[0].trim(0,st[1].stats.endtime - ((len_diff)*st[0].stats.delta))
            print('East component is {} points longer than North component, trimming'.format(len_diff))
        elif  st[1].data.size  >  st[0].data.size : # st[1] longer so trim from its end
            st[1].trim(0,st[1].stats.endtime - ((len_diff)*st[1].stats.delta))
            print('North component is {} points longer than East component, trimming'.format(len_diff))

    # print(st[0].data.size,st[1].data.size)

    if st[0].data.size %2 == 0 or st[1].data.size %2 == 0 : # tests to see if there is an even nukber of points
        # print("Even number of points, trimming by 3")
        st = st.trim(0,st[0].stats.endtime - 3*(st[0].stats.delta))
        # print(st[0].data.size,st[1].data.size)

    st.filter("bandpass",freqmin= f_min, freqmax= f_max,corners=2,zerophase=True) # Zerophase bandpass filter of the streams
    pair = sw.Pair(st[1].data,st[0].data,delta = st[0].stats.delta)
    return pair

def eigen_plot(eign,fig,**kwargs):
    """
    Function broadly copied from SplitwavePy for plotting splitting measurements from an eigm file.
    The difference being that this verison also returns the figure handle so the figure can be manipulated to add an interactive quality/re-do picker
    """
    # setup figure and subplots
    #fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(2, 3,width_ratios=[1,1,2])
    ax0 = plt.subplot(gs[0,0])
    ax1 = plt.subplot(gs[0,1])
    ax2 = plt.subplot(gs[1,0])
    ax3 = plt.subplot(gs[1,1])
    ax4 = plt.subplot(gs[:,2])

    # data to plot
    d1 = eign.data.chop()
    d1f = eign.srcpoldata().chop()
    d2 = eign.data_corr().chop()
    d2s = eign.srcpoldata_corr().chop()

    # flip polarity of slow wave in panel one if opposite to fast
    # d1f.y = d1f.y * np.sign(np.tan(eign.srcpol()-eign.fast))

    # get axis scaling
    lim = np.abs(d2s.data()).max() * 1.1
    ylim = [-lim,lim]

    # original
    d1f._ptr(ax0,ylim=ylim,**kwargs)
    d1._ppm(ax1,lims=ylim,**kwargs)
    # corrected
    d2s._ptr(ax2,ylim=ylim,**kwargs)
    d2._ppm(ax3,lims=ylim,**kwargs)

    # error surface
    if 'vals' not in kwargs:
        # kwargs['vals'] = (eign.lam1 - eign.lam2) / eign.lam2
        # kwargs['title'] = r'$(\lambda_1 - \lambda_2) / \lambda_2$'
        kwargs['vals'] = eign.lam1 / eign.lam2
        kwargs['title'] = r'$\lambda_1 / \lambda_2$'

    # add marker and info box by default
    if 'marker' not in kwargs: kwargs['marker'] = True
    if 'info' not in kwargs: kwargs['info'] = True
    if 'conf95' not in kwargs: kwargs['conf95'] = True
    eign._psurf(ax4,**kwargs)
    # neaten
    plt.tight_layout()


def split_measure(pair,SKS):
    """
    Function for Picking the window for a provded pair object and then measure the splitting
    """
    pair.plot(pick=True,marker = SKS) # Plots the window picker.
    split = sw.EigenM(pair,lags=(4,) )

    fig = plt.figure(figsize=(12,6))
    eigen_plot(split,fig) #Call function eigen_plot to plot splititng messuremnt

    cid = fig.canvas.mpl_connect('key_press_event',interact)
    plt.show(fig)
    fig.canvas.mpl_disconnect(cid)
    return split, pair.wbeg(), pair.wend(),fig

def interact(event):
    """
    Function for quality assignment.
    Responses differe depending on kep press

    """
    if event.key in ['a','b','c']:
        quality.append(event.key)
        plt.close()
    elif event.key is 'n':
        quality.append('n')
        plt.close()
    elif event.key is 'x':
        quality.append('x')
        plt.close()
    elif event.key == 'r':
        plt.close()
        split_measure(pair_glob)
    else:
        print('Invalid key_press_event, please press a,b,c,n,r or x')

# def null():
#     """
#     Functon to confirm the observation of a null event.
#     """
#     var = input("Null Measurement (y/n)")
#     if var is 'y':
#         quality.append('null')
#         print("Null Confirmed")
#         plt.close()
#     elif var is 'n':
#         print("Ok, repeating measurement")
#         split_measure(pair_glob)
#     else:
#         print("Invalid response, try again")
#         null()

def Jacks_SKS_RAW(station):
    """
    Function to read in Jack Walpoles Raw data for comparison for a given station
    """
    raw = pd.read_csv("../Data/Jacks_SKS_RAW.txt",delim_whitespace=True)
    JACK = raw[(raw['STAT'] == station) ]
    JACK = JACK.reset_index()
    del JACK['index']

    return JACK

def split_match(date,time,station):
    """
    Function to find and extract the measuremnt from Jack Walpoles splititng data for the same event.
    This matching is done by finding an entry with the same date stamp. Inital testing has shown this to be a unique identifier.
    station MUST be a string of a station code
    date MUST be a int/float of with the format yyyyjjj where j is julian day
    """
#   -------
#   First we need to read in the splitting observations made by Jack Walpole.
#   We also slice out splitting observations just from the station of interest and then reset the indexing so WL_split's indicies start from [0]
#   -------
    WL_split = Jacks_SKS_RAW(station)
#   -------
#   Using a Pandas DataFrame we can slice out any rows that match out dare stamp
#   The the iloc function is used to extract the requisite values (Here this is trivial as match should be a single row dataframe, but the values still need to be extracted this way)
#
    match = WL_split[(WL_split['DATE'] == date)] # slicies rows in WL_split that have the same datestamp as date. In theory this should return a single row DataFrame
    if len(match) == 1:
        (fast,dfast,tlag,dtlag,wbeg,wend) = (match.iloc[0]['FAST'],match.iloc[0]['DFAST'],match.iloc[0]['TLAG'],match.iloc[0]['DTLAG'],match.iloc[0]['WBEG'],match.iloc[0]['WEND'])

    elif len(match) == 0:

        print("The provided datestamp {} does not match any obervations made by JW".format(date))
        (fast,dfast,tlag,dtlag,wbeg,wend) = ('NaN','NaN','NaN','NaN','40','80')
    else:

        print("There has been more than one match, now testing by timestamp also!\n")
        time_test = int(str(time).zfill(6)[0:4]) #Jacks timestamps are only hhmm so I need to strip off the seconds from my timestamps. WARNING it is possible my timestamps are different to Jacks!!
        print('My timestamp {}, Jacks timestamp {}'.format(time_test,match.iloc[0]['TIME']))
        match2 = WL_split[(WL_split['DATE'] == date) & (WL_split['TIME'] == time_test)]
        # print(match2)
        (fast,dfast,tlag,dtlag,wbeg,wend) = (match.iloc[0]['FAST'],match.iloc[0]['DFAST'],match.iloc[0]['TLAG'],match.iloc[0]['DTLAG'],match.iloc[0]['WBEG'],match.iloc[0]['WEND'])

        if len(match2) == 0: #If there is still no match
            (fast,dfast,tlag,dtlag,wbeg,wend) = ('NaN','NaN','NaN','NaN','NaN','NaN')
            print("No match found")
    return fast,dfast,tlag,dtlag,wbeg,wend

def diag_plot(file,title1):
    """
    Function to make diagnostice plots for a given file of splitting measuremtns
    """
    data = pd.read_csv(file,delim_whitespace=True)
    a = data['FAST']
    d = data.index[np.isnan(data['FAST']) == True].tolist() # Find any rows which contain NaNs
    data = data.drop(d)
    data = data[(data.QUAL != 'x')]
    fig,axs = plt.subplots(2, 2,sharex='col',figsize=(10,10))

    plt.subplot(221)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].FAST,yerr=data[(data.QUAL == 'n')].DFAST,fmt='kx',elinewidth=0.5,label='Null')
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].FAST,yerr=data[(data.QUAL != 'n')].DFAST,fmt='ko',elinewidth=0.5,label='Split')
    plt.legend(loc=2)

    plt.ylabel('Fast Direction (deg)')
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title('{} - Fast Direction'.format(title1))

    plt.subplot(223)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_FAST,yerr=data[(data.QUAL == 'n')].WL_DFAST,fmt='kx',elinewidth=0.5)
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_FAST,yerr=data[(data.QUAL != 'n')].WL_DFAST,fmt='ko',elinewidth=0.5)
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title('Jacks(Sheba) - Fast Direction')
    plt.xlabel('Back Azimuth')
    plt.ylabel('Fast Direction (deg)')

    plt.subplot(222)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].TLAG,yerr=data[(data.QUAL == 'n')].DTLAG,fmt='kx',elinewidth=0.5)
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].TLAG,yerr=data[(data.QUAL != 'n')].DTLAG,fmt='ko',elinewidth=0.5)
    plt.ylabel('Tlag (s)')
    plt.ylim([0,4])
    plt.title('{} - Lag Time'.format(title1))

    plt.subplot(224)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_TLAG,yerr=data[(data.QUAL == 'n')].WL_DTLAG,fmt='kx',elinewidth=0.5)
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_TLAG,yerr=data[(data.QUAL != 'n')].WL_DTLAG,fmt='ko',elinewidth=0.5)
    plt.ylim([0,4])
    plt.ylabel('Tlag (s)')
    plt.xlabel('Back Azimuth')
    plt.title('Jacks(Sheba) - Lag Time')


    plt.tight_layout()
    plt.show()
