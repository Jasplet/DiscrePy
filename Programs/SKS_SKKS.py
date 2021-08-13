#!/usr/bin/env python

"""
Script to measure SKS/SKKS splitting on SAC data
"""

import glob
import numpy as np
import splitwavepy
import obspy
from obspy import UTCDateTime
from obspy.taup import TauPyModel
model = TauPyModel('iasp91')

# josephs results file
jo_results = np.genfromtxt('/Users/glyjw/Joseph/SKS_SKKS/SKS_SKKS_110/SKS_SKKS_110_pairs.sdb', names=True, dtype=None)

# the directory holding the data
directory = '/Users/ja17375/Shear_Wave_Splitting/Sheba/

# place to save measurements
outputdir = '/Users/glyjw/Joseph/SKS_SKKS/SKS_SKKS_110/splitting/'

# the main function
def measure_sks_skks(filestem):

    st = obspy.read(directory + filestem + '[N,E].sac')
    # apply a bandpass filter
    st.filter("bandpass",freqmin=0.01,freqmax=0.5)

    # estimate SKS and SKKS arrival times
    info = st[0].stats.sac
    evlat, evlon = info['evla'], info['evlo']
    stlat, stlon = info['stla'], info['stlo']
    evdepth = info['evdp']

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
    splitstem = filestem.split('_')
    line = jo_results[ (jo_results['STAT'] == splitstem[0].encode('ascii')) &
                (jo_results['DATE'] == int(splitstem[1])) ]
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
    data.EigenM(lags=(4,), bootstrap=True).save(outputdir + filestem + '.sks.eigm')
    data.TransM(lags=(4,), pol=baz, bootstrap=True).save(outputdir + filestem + '.sks.trnm')
    data.XcorrM(lags=(4,), bootstrap=True).save(outputdir + filestem + '.sks.xcrm')

    # measure skks splitting
    data.set_window(30+(wbeg_skks-wbeg_sks), 30+(wend_skks-wbeg_sks))
    data.EigenM(lags=(4,), bootstrap=True).save(outputdir + filestem + '.skks.eigm')
    data.TransM(lags=(4,), pol=baz, bootstrap=True).save(outputdir + filestem + '.skks.trnm')
    data.XcorrM(lags=(4,), bootstrap=True).save(outputdir + filestem + '.skks.xcrm')

# Loop through files
for f in glob.glob(directory + '/*BHZ.sac'):
    filestem = f.split('/')[-1][:-5]
    print('processing ' + filestem)
    measure_sks_skks(filestem)
