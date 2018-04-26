#! /usr/bin/env python
#  Inspecter.py
#  A program to loop through and visiualy inspect/ assign data qauility and qulaitvative matching to SKS-SKKS pairs
#### Import Modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import cartopy.crs as cart
#import cartopy
import matplotlib.gridspec as gridspec
import obspy
from obspy import taup
import os
import sys


def measure_sks_skks(filestem,outputdir,):
    """
    Function to measure SWS using the eigenvalue, transverse minimisation and cross correlation methods
    This is done using the devel release of splitwave py (Needs to run from the splitwavepy environment)
    """
    st = obspy.read(directory + filestem + '[N,E].sac')
    f_check_SKS = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/{}/{}/SKS/{}SKS_sheba.final_result'.format(run,st[0].stats.station,phase,filestem)
    f_check_SKKS =  '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/{}/{}/SKKS/{}SKKS_sheba.final_result'.format(run,st[0].stats.station,phase,filestem)
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

    wbeg_sks = evtime + line['WBEG_SKS']
    wend_sks = evtime + line['WEND_SKS']
    wbeg_skks = evtime + line['WBEG_SKKS']
    wend_skks = evtime + line['WEND_SKKS']
    st.trim(wbeg_sks - 30, wend_skks + 30)

    # put data into splitwavepy object
    data = splitwavepy.Data(st[1].data, st[0].data, delta=st[0].stats.delta)
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
    data.EigenM(lags=(4,), bootstrap=True).save('{}/{}.sks.eigm'.format(outputdir,filestem))
    data.TransM(lags=(4,), pol=baz, bootstrap=True).save('{}/{}.sks.trnm'.format(outputdir,filestem))
    data.XcorrM(lags=(4,), bootstrap=True).save('{}/{}.sks.xcrm'.format(outputdir,filestem))

    # measure skks splitting
    data.set_window(30+(wbeg_skks-wbeg_sks), 30+(wend_skks-wbeg_sks))
    data.EigenM(lags=(4,), bootstrap=True).save('{}/{}.skks.eigm'.format(outputdir,filestem))
    data.TransM(lags=(4,), pol=baz, bootstrap=True).save('{}/{}.skks.trnm'.format(outputdir,filestem))
    data.XcorrM(lags=(4,), bootstrap=True).save('{}/{}.skks.xcrm'.format(outputdir,filestem))

# Loop through files

if __name__ == '__main__':
    print('Hello, this is inspecter.py')

    # Read in the matched SKS SKKS pairs
    pairs = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}'.format(sys.argv[1]))
    # Define which directory the SAC files are saved in
    data_dir = '/Users/ja17375/Shear_Wave_Splitting/Sheba'
    # Define Output dir
    outputdir = 'Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/QA'
    # Test to see if the out dir exists (and if not make it!)
    if os.path.isdir(outputdir) is False:
        os.makedirs(outputdir)


    evt_list = sys.argv[2]
    file_list ='/Users/ja17375/Shear_Wave_Splitting/Data/{}'.format(evt_list)
    # echo out where I expect the staiton list to be
    print('Processing Data from the Downloaded Event List {}'.format(file_list))
    files = []
    with open(file_list,'r') as reader:
        for line in reader.readlines():
            f = line.strip('\n')
            files.append(f)

    for file in files:
        filestem = file[55:].strip('/')

        splitstem = filestem.split('_')
        line = pairs[(pairs['STAT'] == splitstem[0].encode('ascii')) &
                    (pairs['DATE'] == int(splitstem[1])) ]

        measure_sks_skks(filestem)
