#! /anaconda3/envs/splitwavepy/bin/python

# Benchmark swp measurements against sheba results
### Imports
import numpy as np
import pandas as pd
import obspy as ob
import splitwavepy as sw
import matplotlib.pyplot as plt
import glob
import Split_Measure as sm
# Load sheba results

sheba = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Jacks_Split_SKS_sheba_results.sdb',delim_whitespace=True)

# Right now we only care about SKS results so lets parse them to new variables
lag_sheba,dlag_sheba = sheba.TLAG.values, sheba.DTLAG.values
fast_sheba, dfast_sheba = sheba.FAST.values, sheba.DFAST.values
wbeg,wend = sheba.WBEG.values, sheba.WEND.values
# We also need to back-azimuth vlaue to plot against
baz = sheba.BAZ.values
# Read the SplitWavePy results
files = glob.glob('/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/SplitWavePy/*sks.eigm')

lag_swp,dlag_swp = [],[]
fast_swp,dfast_swp = [],[]

for i,f in enumerate(files):
    # First extract filestem (E.G. NEW_2000100_0001_ )
    fstem = f.split('/')[-1].strip('_sks.eigm')
    stat = fstem.split('_')[0] # Isolates station code as we need that seperatly
    try:
        st = sm.read_sac('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}*'.format(stat,fstem))

    except Exception:
        print(Exception)
        fstem2 = fstem[0:-1]
        st = sm.read_sac('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}*'.format(stat,fstem2))

    st.filter('bandpass',freqmin=0.01,freqmax=0.5)

    attribs = st.stats.sac

    times = sm.model_traveltimes(st[0],'SKS') # Returns traveltime in UTC, origin time (UTC) and relative arrivals time

    #Set Windows and trim. TImes[1] should be event origin time
    start = wbeg[i] + times[1]
    end = wend[i] + times[1]
    st.trim(start-30, end+30)
    
    # Create new SWP data object
    sw.Data(st.select(channel='BHN')[0].data, st.select(channel='BHE')[0].data, delta=st.select(channel='BHN')[0].stats.delta)
    # add the SAC header info
    data.sacinfo = st[0].stats.sac)
