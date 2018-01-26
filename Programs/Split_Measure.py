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

def splitting(station,switch,files):
    """
    Measures SKS splitting for all streams listed in a ttext file at the provided path. These streams must be saved as SAC files.abs
    This function is the primary part of this module/package/script/thing, all the pther functions support this one.

    inpath - string contaiing the path of the textfile which contains the sac file to be read

    switch - optional kwarg to specify if you want to manually window the data or use a set of windows (Walpoles windows are availbale for use by default)
    """
    outfile = output_init(station)
    if switch == 'off':
        ui = input('Please enter file name (minus extension) for the eigm files: \n')
        filename = '/Users/ja17375/Python/SKS_Splitting/Eigm_Files/{}'.format(ui)
        # plot_on = input('Do you want to plot the eigm surfaces? (y/n) \n')
    else:
        pass
    with open(files,'r') as reader: # NEW_read_stream.txt is a textfile containing filenames of streams which have been read and saved by Split_Read for this station. s
        for line in reader.readlines():
            line.strip('\n')
            st_id = line[0:-7]
            st = read_sac(st_id+'*')
            # Intialise some global variables which I need to pass things between fucntions (this is probably not be best way to do things but it works!)
            global pair_glob
            global quality
            if st != False: # i.e. if the stream is sufficiently populated and has been read.
#               if eigm != False (if there is already an eigm for this event)
                #   a = input('There is already the measuremnt {} for this event, do you want to remseasure? (y/n)'.format(eigm))

############### Measuring Start Here! #####################
                SKS_UTC, t0, SKS = model_SKS(st[0]) # Returns SKS arrival as a UTCDateTime object, event origin time and SKS arrival relative to the origin time
                # print(t0)
                # print('SKS_UTC ={}'.format(SKS_UTC))
                quality = [] # variable to hold Callback key entries for estimated quality of splitting measurements
                date,time = int(str(t0.year)+str(t0.julday).zfill(3)),int(str(t0.hour).zfill(2)+str(t0.minute).zfill(2)+str(t0.second).zfill(2)) #creates time and date stamps
                pair = st_prep(st = st, f_min = 0.01,f_max = 0.5)
                # print(st[0].stats.starttime)
                pair_glob = pair
                # print(' Predicted SKS arrival is at {:5.3f}\n'.format(rel_SKS))

    #           --------------
    #           Now lets find what splitting Jack Walpole measured for this event
    #           --------------
                (wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend) = split_match(date,time,station)

                if switch == 'on': # If manual windowing is on
                    split, wbeg, wend,fig = split_measure(pair,SKS)
                    eig_file = '{}{}_{:07d}_{:06d}.eigm'.format('/Users/ja17375/Python/SKS_Splitting/Eigm_Files/',station,date,time)
                    plt.savefig('{}{}_{:07d}_{:06d}'.format('/Users/ja17375/Python/SKS_Splitting/Figures/Eigm_Surface/',station,date,time))
                    plt.close()

                elif switch == 'off': #Manual windowing is off. For now this will just mean Jacks windows will be used. Eventually add automation or support for entering windows.
                    wbeg,wend = wl_wbeg,wl_wend

                    pair.set_window(start=wl_wbeg,end=wl_wend) # Sets window to that of Jacks
                    split = sw.EigenM(pair,lags=(4,) ) # Measures splitting

                    fig = plt.figure(figsize=(12,6))
                    eigen_plot(split,fig)
                    plt.savefig('{}/{}_{:07d}_{:06d}'.format('/Users/ja17375/Python/SKS_Splitting/Figures/Eigm_Surface',ui,date,time))
                    plt.close()
                    quality = 'w'
                    eig_file = '{}_{:07d}_{:06d}.eigm'.format(filename,date,time)

                # print(t0+wbeg,t0+wend)
                split.save(eig_file) # Saves splitting measurements
                # if quality is not ('x'): #If the quality attribute is not bad (indicated by x)
                print('WBEG: {}, WEND: {}'.format(wbeg,wend))

                attrib = ['stla','stlo','evla','evlo','evdp','gcarc','baz'] #SAC attribute values that I want to extract and print later
                meas = [wbeg, wend, split.fast, split.dfast, split.lag, split.dlag,wl_fast,wl_dfast,wl_tlag,wl_dtlag,wl_wbeg,wl_wend ]
                stats = [st[0].stats.sac[i] for i in attrib]
                org = [station,date,time]

                outfile.write('{} {:07d} {:06d} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:05.2f} {:06.2f} {:06.2f} {:4.2f} {:4.2f} {:4.2f} {:4.2f} {:5.3f} {:4.2f} {} {} {} {} {} {} {}\n'.format(*org,*stats,*meas,quality[0]))
                save_sac(st,quality[0],date,time,wbeg,wend,switch)
                # else:
                #     meas = ['NaN','NaN','NaN','NaN','NaN','NaN']
                #     stats = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN']
                #     org = ['NEW',date, time]

                    # outfile.write('{} {:07.0d} {:06.0d} {} {} {} {} {} {} {} {} {} {} {} {} {} {} \n'.format(*org,*stats,*meas,quality[0])
            else:
                meas = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN']
                stats = ['NaN','NaN','NaN','NaN','NaN','NaN','NaN']
                org = [station,'NaN','NaN']
                quality = ['NaN']
                print('No stream for event',line[0:-7])
                outfile.write('{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16} {17} {18} {19} {20} {21} {22}\n'.format(*org,*stats,*meas,quality[0]))


    outfile.close()

def output_init(station):
    """
    Initialises output variables and output textfile

    station - string containing the station code
    """
    default_out = '/Users/ja17375/Python/SKS_Splitting/Measurements/{}_Splitting.txt'.format(station) #Default output filename
    #default_out = '/Users/ja17375/Python/SKS_Splitting/Measurements/{}_Splitting_JW_Windows.txt'.format(station)
    if os.path.isfile(default_out):
        #Default file exists! Request user permission to overwrite
        ovr = user_in('c1',station)
        if ovr == 'y':
            outfile = open(default_out,'w+')

        elif ovr == 'n':
            new_out = user_in('c2',station)
            outfile = open('/Users/ja17375/Python/SKS_Splitting/Measurements/{}.txt'.format(new_out),'w+')
    elif os.path.isfile(default_out) == False:
        print('{} does not exist and will be created'.format(default_out))
        outfile = open(default_out,'w+')

    outfile.write('STAT DATE TIME STLA STLO EVLA EVLO EVDP GCARC BAZ WBEG WEND FAST DFAST TLAG DTLAG WL_FAST WL_DFAST WL_TLAG WL_DTLAG WL_WBEG WL_WEND QUAL\n')
    return outfile

def user_in(case,station):
    if case == 'c1':
        a = input('The file {}_Splitting.txt already exists, would you like to overwrite? (y/n)\n'.format(station))
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
    try:
        st = obspy.core.stream.read(st_id) # Reads all traces that match st_id
        if len(st) == 3: #Check is there are 3 traces in the stream (East,North and Vertical)
            return st
        else:
            return False
    except Exception:
        return False


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
            path = '/Users/ja17375/Python/SKS_Splitting/Data/Proccessed_Streams'
            tr2 = tr.trim(t0 + wbeg, t0 + wend)
            tr2.write('{}/{}_{}_{:07d}_{:06d}_{}.sac'.format(path,stat,qual,date,time,ch), format='SAC')
    elif switch == 'off':
        pass

def model_SKS(tr):
    """
    Function to run TauP traveltime models for the SKS phase.
    Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object

    tr - trace object for which SKS arrival time will be predicted
    """
    model = obspy.taup.tau.TauPyModel(model="iasp91")
    SKS = model.get_travel_times(tr.stats.sac.evdp,tr.stats.sac.gcarc,["SKS"])[0].time
    t0 = tr.stats.starttime # Start time of stream for event. This should be the event origin time.
    SKS_UTC = obspy.core.utcdatetime.UTCDateTime(t0 + SKS)# SKS arrival time relative to trace start as a UTCDateTime object

    return SKS_UTC, t0, SKS


def st_prep(st,f_min,f_max):
    """
    Prepares Stream for spltting analysis (by bandpass filtering and trimming) and then creates the Pair object
    """
    if st[0].stats.npts %2 == 0: # tests to see if there is an even nukber of points
        st = st.trim(st[0].stats.starttime,st[0].stats.endtime - st[0].stats.delta)

    st.filter("bandpass",freqmin= f_min, freqmax= f_max,corners=2,zerophase=True) # Zerophase bandpass filter of the streams
    # st.trim(starttime = SKS - trim, endtime = SKS + trim)
    # rel_SKS = SKS - st[0].stats.starttime
    return sw.Pair(st[1].data,st[0].data,delta = st[0].stats.delta)

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
    plt.show()

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

def null():
    """
    Functon to confirm the observation of a null event.
    """
    var = input("Null Measurement (y/n)")
    if var is 'y':
        quality.append('null')
        print("Null Confirmed")
        plt.close()
    elif var is 'n':
        print("Ok, repeating measurement")
        split_measure(pair_glob)
    else:
        print("Invalid response, try again")
        null()

def interact(event):
    """
    Function for quality assignment.
    Responses differe depending on kep press

    """
    if event.key in ['a','b','c']:
        quality.append(event.key)
        plt.close()
    elif event.key is 'n':
        null()
    elif event.key is 'x':
        quality.append('x')
        plt.close()
    elif event.key == 'r':
        plt.close()
        split_measure(pair_glob)
    else:
        print('Invalid key_press_event, please press a,b,c,n,r or x')

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
    fig,axs = plt.subplots(2, 2,sharex='col',figsize=(10,10))

    plt.subplot(221)
    plt.errorbar(data['BAZ'],data['FAST'],yerr=data['DFAST'],fmt='o',elinewidth=0.5)

    plt.ylabel('Fast Direction (deg)')
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title('{} - Fast Direction'.format(title1))

    plt.subplot(223)
    plt.errorbar(data['BAZ'],data['WL_FAST'],yerr=data['WL_DFAST'],fmt='ro',elinewidth=0.5)
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title('Jacks(Sheba) - Fast Direction')
    plt.xlabel('Back Azimuth')
    plt.ylabel('Fast Direction (deg)')

    plt.subplot(222)
    plt.errorbar(data['BAZ'],data['TLAG'],yerr=data['DTLAG'],fmt='o',elinewidth=0.5)
    plt.ylabel('Tlag (s)')
    plt.ylim([0,4])
    plt.title('{} - Lag Time'.format(title1))

    plt.subplot(224)
    plt.errorbar(data['BAZ'],data['WL_TLAG'],yerr=data['WL_DTLAG'],fmt='ro',elinewidth=0.5)
    plt.ylim([0,4])
    plt.ylabel('Tlag (s)')
    plt.xlabel('Back Azimuth')
    plt.title('Jacks(Sheba) - Lag Time')

    plt.tight_layout()
    plt.show()
