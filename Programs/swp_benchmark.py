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
from scipy.stats import gaussian_kde

###################################

def benchmark(path):
    'Make swp measurements to correspond with those from sheba '
    # Load sheba results
    sheba = pd.read_csv(path,delim_whitespace=True)
    # Right now we only care about SKS results so lets parse them to new variables
    if path.split('.')[1] == 'sdb':
        sheba_results = sheba.FAST.values, sheba.DFAST.values,sheba.TLAG.values, sheba.DTLAG.values
        wbeg,wend = sheba.WBEG.values, sheba.WEND.values

    elif path.split('.')[1] == 'pairs':
        sheba_results = sheba.FAST_SKS.values, sheba.DFAST_SKS.values,sheba.TLAG_SKS.values, sheba.DTLAG_SKS.values
        wbeg,wend = sheba.WBEG_SKS.values, sheba.WEND_SKS.values

    else:
        print('Unsupported File Type')
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
            st = ob.read('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}_*_BH[N,E].sac'.format(stat,fstem2))
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

    fig,axs = plt.subplots(2, 2,sharex='col',figsize=(10,10))

    plt.subplot(221)
    #plt.errorbar(baz,sheba[0],yerr=sheba[1],fmt='kx',elinewidth=0.5,label='Sheba')
    plt.plot(baz,sheba[0],'kx')
    #plt.legend(loc=2)

    plt.ylabel('Fast Direction (deg)')
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title(r'$\phi$ from Sheba')

    plt.subplot(222)
    #plt.errorbar(baz,sheba[2],yerr=sheba[3],fmt='kx',elinewidth=0.5,label='Sheba')
    plt.plot(baz,sheba[2],'kx')

    plt.ylim([0,4])
    plt.yticks(np.arange(0,4,0.5))
    plt.title(r'$\delta t$ from Sheba')
    plt.xlabel('Back Azimuth (deg)')
    plt.ylabel(r'$\delta t$')

    plt.subplot(223)
    #plt.errorbar(baz,swp[0],yerr=swp[1],fmt='bx',elinewidth=0.5,label='SplitWavePy')
    plt.plot(baz,swp[0],'bx')

    #plt.legend(loc=2)
    plt.ylabel('Fast Direction (deg)')
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title(r'$\phi$ from SplitwavePy')

    plt.subplot(224)
    #plt.errorbar(baz,swp[2],yerr=swp[3],fmt='bx',elinewidth=0.5,label='SplitWavePy')
    plt.plot(baz,swp[2],'bx')

    plt.title(r'$\delta t$ from SplitwavePy')
    plt.xlabel('Back Azimuth (deg)')
    plt.ylabel(r'$\delta t$')
    plt.ylim([0,4])
    plt.tight_layout()
    plt.show()

def splitting_kde(lag ,fast):
    """
    Function to perform kernal density estimation for set of splitting results.
    """

    L,F = np.mgrid[ 0:4:500j,-90:90:500j]
    positions = np.vstack([L.ravel(),F.ravel()])
    splitting = np.vstack([lag,fast])
    kernal = gaussian_kde(splitting)
    S = kernal.evaluate(positions)
    surf = S.reshape(L.shape)
    # When the surface is returend we need to transpose it so that it is correctly oriented (for some reason)
    return np.transpose(surf)

def sheba_v_swp(sheba,swp):
    '''
    Plot sheba fast and lag results against those from splitwavepy
    '''
    fig = plt.subplots(nrows=2,ncols=1,figsize=(9,9))
    plt.subplot(211)

    # Plot sheba v swp FAST direction
    plt.errorbar(sheba[0],swp[0],xerr=sheba[1],yerr=swp[1],fmt='k.',label='Results',elinewidth=0.5)
    plt.plot(np.linspace(-90,90,50),np.linspace(-90,90,50),'k-.',label='Expected')
    plt.legend(loc=0)
    plt.title(r'Comparison of $\phi$ for SHEBA and SplitWavePy')
    plt.xlabel('Results from sheba')
    plt.ylabel('Results from SplitWavePy')
    plt.yticks(np.arange(-90,91,30))
    plt.xticks(np.arange(-90,91,30))

    # Plot sheba v swp lag times
    plt.subplot(212)

    plt.errorbar(sheba[2],swp[2],xerr=sheba[3],yerr=swp[3],fmt='k.',label='Results',elinewidth=0.5)
    plt.plot(np.linspace(0,4,50),np.linspace(0,4,50),'k-.',label='Expected')
    plt.legend(loc=0)
    plt.title(r'Comparison of $\delta t$ for SHEBA and SplitWavePy')
    plt.xlabel(r'$\delta t$ from sheba')
    plt.ylabel(r'$\delta t$ from SplitWavePy ')

    plt.show()

def fast_lag_plot(sheba_lag,sheba_fast,swp_lag,swp_fast):
    """
    Plot fast v lag for swp and sheba
    """


    fig,axs = plt.subplots(2,1,sharex='col',figsize=(8,8))
    # Plot tlag on x and fast on y ( so it looks a bit like a lambda surface)
    plt.subplot(211)
    # For sheba results
    f = np.linspace(-90,90,500)
    l = np.linspace(0,4,500)
    surf = splitting_kde(sheba_lag,sheba_fast)
    #plt.imshow(surf,extent=[0,4,-90,90],aspect='auto')

    plt.contour(l,f,surf)
    plt.plot(sheba_lag,sheba_fast,'k.',markersize=1.5)
    plt.xlim([0,4])
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title(r'$\phi$ v $\delta t$ for SHEBA results')

    plt.subplot(212)
    # For splitwave py results
    surf = splitting_kde(swp_lag,swp_fast)
    plt.contour(l,f,surf)
    plt.plot(swp_lag,swp_fast,'k.',markersize=1.5)
    plt.xlim([0,4])
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title(r'$\phi$ v $\delta t$ for SplitWavePy results')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':

    sheba,swp,baz = benchmark('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Jacks_Split_SKS_sheba_results.sdb')
    diag_plot(sheba,swp,baz)
    fast_lag_plot(sheba[2],sheba[0],swp[2],swp[0])
    # sheba[0] should hold fast results from sheba
    # sheba[2] should hold lag results from sheba
    # The same holds for swp
