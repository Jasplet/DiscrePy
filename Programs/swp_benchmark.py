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
from datetime import datetime
###################################

def benchmark(path):
    'Make swp measurements to correspond with those from sheba '
    # Load sheba results
    sheba = pd.read_csv(path,delim_whitespace=True)
    # Right now we only care about SKS results so lets parse them to new variables
    sheba_results = sheba.FAST.values, sheba.DFAST.values,sheba.TLAG.values, sheba.DTLAG.values

    wbeg,wend = sheba.WBEG.values, sheba.WEND.values
    # We also need to back-azimuth vlaue to plot against
    baz = sheba.BAZ.values

    # Read the SplitWavePy results
    files = glob.glob('/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/SplitWavePy/*sks.eigm')

    lag_swp,dlag_swp = [],[]
    fast_swp,dfast_swp = [],[]

    for i,f in enumerate(sheba.DATE.values):
        print('It {}, time is {} '.format(i,str(datetime.now())))
        # First get the right DATE,TIME and STATION from sehab results
        date,time,stat = sheba.DATE[i], sheba.TIME[i], sheba.STAT[i]
        fstem = '{}_{}_{}'.format(stat,date,time)
        try:
            st = ob.read('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}??_BH[N,E].sac'.format(stat,fstem))

        except Exception:
            print(Exception)
            fstem2 =  '{}_{}'.format(stat,date)
            st = ob.read('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}_??????_BH[N,E].sac'.format(stat,fstem2))
        #print(st)
        st.filter('bandpass',freqmin=0.01,freqmax=0.5)

        attribs = st[0].stats.sac

        times = sm.model_traveltimes(st[0],'SKS') # Returns traveltime in UTC, origin time (UTC) and relative arrivals time

        #Set Windows and trim. TImes[1] should be event origin time
        start = times[1] + wbeg[i]
        end = times[1] + wend[i]
        st.trim(start-30, end+30)

        # Create new SWP data object
        data = sw.Data(st.select(channel='BHN')[0].data, st.select(channel='BHE')[0].data, delta=st.select(channel='BHN')[0].stats.delta)
        # add the SAC header info
        data.sacinfo = st[0].stats.sac

        data.set_window(30,30+(end-start))
        swp = data.SC(lags=(4,))
        lag_swp.append(swp.lag)
        dlag_swp.append(swp.dlag)
        fast_swp.append(swp.fast)
        dfast_swp.append(swp.dfast)

    swp_results = fast_swp,dfast_swp,lag_swp,dlag_swp
    print('{} sheba results, {} swp results'.format(len(sheba_results[0]),len(swp_results[0])))
    return sheba_results,swp_results,baz

def diag_plot(sheba,swp,baz):
    """
    Function to make diagnostice plots for a given file of splitting measuremtns
    """

    fig,axs = plt.subplots(2, 1,sharex='col',figsize=(10,10))

    plt.subplot(211)
    plt.errorbar(baz,sheba[0],yerr=sheba[1],fmt='kx',elinewidth=0.5,label='Sheba')
    plt.errorbar(baz,swp[0],yerr=swp[1],elinewidth=0.5,label='SplitWavePy')
    plt.legend(loc=2)

    plt.ylabel('Fast Direction (deg)')
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title(r'Comparison of $\phi$)')

    plt.subplot(212)
    plt.errorbar(baz,sheba[2],yerr=sheba[3],fmt='kx',elinewidth=0.5,label='Sheba')
    plt.errorbar(baz,swp[2],yerr=swp[3],elinewidth=0.5,label='SplitWavePy')
    plt.ylim([0,4])
    plt.yticks(np.arange(0,4,0.5))
    plt.title(r'Comparison of $\delta t$')
    plt.xlabel('Back Azimuth (deg)')
    plt.ylabel(r'$\delta t$')


    plt.tight_layout()
    plt.show()

if __name__ == '__main__':

    sheba,swp,baz = benchmark('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Jacks_Split_SKS_sheba_results.sdb')
    diag_plot(sheba,swp,baz)
