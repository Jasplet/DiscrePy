#! /usr/bin/env python
#  Inspecter.py
'''A program to loop through and visiualy inspect/ assign data qauility and qulaitvative matching to SKS-SKKS pairs'''
#### Import Modules
import pandas as pd
import numpy as np
import splitwavepy
import obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
import obspy
from obspy import taup
import os
import sys


def measure_sks_skks(file,outputdir,windows):
    """
    Function to measure SWS using the eigenvalue, transverse minimisation and cross correlation methods
    This is done using the devel release of splitwave py (Needs to run from the splitwavepy environment)
    """
    filestem = file[55:]
    try:
        st = obspy.read('{}??_BH[N,E].sac'.format(file))
    except Exception:
        print(Exception)
        st = obspy.read('{}_??????_BH[N,E].sac'.format(file[:-5]))
    f_check_SKS = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/{}/SKS/{}SKS_sheba.final_result'.format(st[0].stats.station,filestem)
    f_check_SKKS =  '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/{}/SKKS/{}SKKS_sheba.final_result'.format(st[0].stats.station,filestem)
    # apply a bandpass filter
    st.filter("bandpass",freqmin=0.01,freqmax=0.5)

    # estimate SKS and SKKS arrival times
    info = st[0].stats.sac
    evlat, evlon = info['evla'], info['evlo']
    stlat, stlon = info['stla'], info['stlo']
    evdepth = info['evdp']

    model = taup.TauPyModel('iasp91')
    arrivals = model.get_travel_times_geo(evdepth,evlat,evlon,stlat,stlon,phase_list=['SKS','SKKS'])
    evtime = UTCDateTime(year=info['nzyear'], julday=info['nzjday'],
                         hour=info['nzhour'], minute=info['nzmin'],
                         second=info['nzsec'], microsecond=1000*info['nzmsec'])
    skstime = evtime + arrivals[0].time
    skkstime = evtime + arrivals[1].time
    delta_t = skkstime - skstime

    # find back azimuth
    baz = info['baz']

    # trim around SKS and SKKS

    wbeg_sks = evtime + windows[0]
    wend_sks = evtime + windows[1]
    wbeg_skks = evtime + windows[2]
    wend_skks = evtime + windows[3]
    st.trim(wbeg_sks - 30, wend_skks + 30)

    # put data into splitwavepy object
    data = splitwavepy.Data(st.select(channel='BHN')[0].data, st.select(channel='BHE')[0].data, delta=st.select(channel='BHN')[0].stats.delta)
    # add the SAC header info
    data.sacinfo = st[0].stats.sac
    #
    # # measure sks splitting
    # data.set_window(30, 30+(wend_sks-wbeg_sks))
    # data.EigenM(lags=(4,)).save(outputdir + filestem + '.sks.eigm')
    # data.TransM(lags=(4,), pol=baz).save(outputdir + filestem + '.sks.trnm')
    # data.XcorrM(lags=(4,)).save(outputdir + filestem + '.sks.xcrm')
    #
    # # measure skks splitting
    # data.set_window(30+(wbeg_skks-wbeg_sks), 30+(wend_skks-wbeg_sks))
    # data.EigenM(lags=(4,)).save(outputdir + filestem + '.skks.eigm')
    # data.TransM(lags=(4,), pol=baz).save(outputdir + filestem + '.skks.trnm')
    # data.XcorrM(lags=(4,)).save(outputdir + filestem + '.skks.xcrm')

    # measure sks splitting
    data.set_window(30, 30+(wend_sks-wbeg_sks))
    data.EigenM(lags=(4,), bootstrap=True).save('{}/{}_sks.eigm'.format(outputdir,filestem))
    data.TransM(lags=(4,), pol=baz, bootstrap=True).save('{}/{}_sks.trnm'.format(outputdir,filestem))
    data.XcorrM(lags=(4,), bootstrap=True).save('{}/{}_sks.xcrm'.format(outputdir,filestem))

    # measure skks splitting
    data.set_window(30+(wbeg_skks-wbeg_sks), 30+(wend_skks-wbeg_sks))
    data.EigenM(lags=(4,), bootstrap=True).save('{}/{}_skks.eigm'.format(outputdir,filestem))
    data.TransM(lags=(4,), pol=baz, bootstrap=True).save('{}/{}_skks.trnm'.format(outputdir,filestem))
    data.XcorrM(lags=(4,), bootstrap=True).save('{}/{}_skks.xcrm'.format(outputdir,filestem))

# Loop through files

if __name__ == '__main__':
    print('Hello, this is SKS_SKKS_qa.py')
    #
    results_dir = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split'
    # Define which directory the SAC files are saved in
    data_dir = '/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files'
    # Define Output dir
    outputdir = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/SplitWavePy'
    # Test to see if the out dir exists (and if not make it!)
    if os.path.isdir(outputdir) is False:
        print('Outdir:{} Exists? {}'.format(outputdir, os.path.isdir(outputdir)))
        os.makedirs(outputdir)

    # Read in the matched SKS SKKS pairs
    pairs = pd.read_csv('{}/{}'.format(results_dir,sys.argv[1]),delim_whitespace=True,converters={'TIME': lambda x: str(x)})

    for i in range(0,len(pairs)):
        file = '{}/{}/{}_{}_{}'.format(data_dir,pairs.STAT[i],pairs.STAT[i],pairs.DATE[i],pairs.TIME[i])
        #print(file[:-5])
        #print('Iteration: {}, File {}'.format(i,file[55:]))
        windows = [pairs.WBEG_SKS[i],pairs.WEND_SKS[i],pairs.WBEG_SKKS[i],pairs.WEND_SKKS[i]]
        measure_sks_skks(file,outputdir,windows)
