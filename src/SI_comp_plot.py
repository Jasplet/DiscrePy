#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 13:32:06 2021

Revisiting some older SKS, SKKS, SYNTH data to make plots for chapter 2 of my thesis
Chapter 2, i.e., global data collection chapter

Wrangles data and makes SI comp. plot, synthetics SNR v splitting params, |BAZ - SPOL| histograms, etc., 

@author: ja17375
"""
import pandas as pd
import matplotlib.pyplot as plt

RES_PATH = '/Users/ja17375/Projects/Corrections'

   
def read_and_cat_synth():
    
    df = pd.DataFrame()
    for spol in ['30', '45', '60']:
        for noise in ['01', '025', '005']:
             sp = pd.read_csv(f'{RES_PATH}/SP{spol}_noise_{noise}_SYNTH_sheba_results.sdb',
                              delim_whitespace=True)
             df = df.append(sp)    
    
    return df

def plot_si(ax, df):
    data = df[df.SNR > 5.0]
    print(len(data))
    s = data[data['Q'] > 0.5]
    n = data[data['Q'] < -0.5]
    rem = data[(data['Q'] >= -0.5) & (data['Q'] <= 0.5)]
    ax.plot(rem['SI(Pa)'],rem['SI(Pr)'],'k.',label='uID') # Plot SI methods for unclear events
    ax.plot(n['SI(Pa)'],n['SI(Pr)'],'b.',label='Null') # Plot SI methods for clear nulls
    ax.plot(s['SI(Pa)'],s['SI(Pr)'],'r.',label='Split') # Plot SI methods for clear splits
    ax.plot([-5,0,5],[-5,0,5],'k--') # add y=x
    return ax

def SI_bowties(sks, skks, scs, synth):
    # Make figure
    fig, axs = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(7,7))
    
    ax1 = axs[0,0]
    plot_si(ax1, sks)
    ax1.set_title('SKS')
    ax1.set_ylabel('SI by projection')

    
    ax2 = axs[0,1]
    plot_si(ax2, skks)
    ax2.set_title('SKKS')
    
    ax3 = axs[1,0]
    plot_si(ax3, scs)
    ax3.set_title('ScS')
    ax3.set_ylabel('SI by projection')
    ax3.set_xlabel('SI by approximation')
    
    ax4 = axs[1,1]
    plot_si(ax4, synth)
    ax4.set_title('Synthetics')
    ax4.set_xlabel('SI by approximation')
    ax4.legend()
    
    ax1.set_xlim([-5,5])
    ax1.set_ylim([-4,4])
    plt.tight_layout()
    plt.savefig('/Users/ja17375/Projects/Corrections/Figs/', dpi=600)

def plot_SI_v_phi_dt(data,save=False):

    fig, axs = plt.subplots(2,2, sharey=True, figsize=(9,9))
    
    # Plot SI(Pa) v tlag on top left axis
    ax1 = axs[0,0]
    C1 = ax1.scatter(data['TLAG'], data['SI(Pa)'], c=data['Q'], cmap='seismic',
                     marker='.')
    ax1.plot([0,2,4],[0,2,4], '--k')
    ax1.plot([0,2,4],[0,-2,-4],'--k')
    # Plot SI(Pr) v tlag on bottom left
    ax2 = axs[1,0]
    ax2.scatter(data['TLAG'], data['SI(Pr)'], c=data['Q'], cmap='seismic',
                     marker='.')
    ax2.plot([0,2,4],[0,2,4], '--k')
    ax2.plot([0,2,4],[0,-2,-4],'--k')
    #Plot SI(Pa) v fast on top right axis
    ax3 = axs[0,1]
    C1 = ax3.scatter(data['FAST'], data['SI(Pa)'], c=data['Q'], cmap='seismic',
                     marker='.')
    fig.colorbar(C1, ax=ax3, label='Q')
    # Plot SI(Pr) v fast on bottom right
    ax4 = axs[1,1]
    C2 = ax4.scatter(data['FAST'], data['SI(Pr)'], c=data['Q'], cmap='seismic',
                     marker='.')
    fig.colorbar(C2, ax=ax4, label='Q')
    ax2.set_xlabel('Delay time (s)')
    ax4.set_xlabel(r'Fast direction $(\degree)$')
    ax1.set_ylabel('SI(Pa)')
    ax2.set_ylabel('SI(Pr)')
    ax1.set_xlim([0,4])
    ax2.set_xlim([0,4])
    ax3.set_xlim([-90,90])
    ax3.set_xticks([-90,-60,-30,0,30,60,90])
    ax4.set_xticks([-90,-60,-30,0,30,60,90])    
    ax4.set_xlim([-90,90])
    if save:
        print('saving')
        plt.savefig('/Users/ja17375/Projects/Corrections/SI_v_splitting_params.png', dpi=600)
    plt.show()

def plot_SI_v_periodfrac():
    ''' Plots plitting intensity versus the dt/T where T is the dominant period of the waveform
       
        This can be only done as synthetics as I know the input dominant frequency = 0.1Hz 
        and thereforet the dominant period is 1/0.1 = 10 s
    '''
    data = read_and_cat_synth()
    frac = data['TLAG'].values / 10 # ratio of dt/T
    
    fig, (ax1, ax2) = plt.subplots(2,1, figsize=(9,7), sharex=True)
    
    C = ax1.scatter(frac, data['SI(Pa)'], c=data['Q'], cmap='seismic',
                marker='.', vmin=-1, vmax=1)
    ax2.scatter(frac, data['SI(Pr)'], c=data['Q'], cmap='seismic',
                marker='.', vmin=-1, vmax=1)
    fig.colorbar(C, ax=[ax1,ax2], label='Q')
    ax1.set_ylabel('SI by Approximation')
    ax2.set_ylabel('SI by Projection')
    ax2.set_xlabel(r'$\delta t$ / T ')

    plt.show()
    
def plot_SI_v_fast(data):
    
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(7,9))
    
    # Plot SI(Pa) v fast on top axis
    C1 = ax1.scatter(data['FAST'], data['SI(Pa)'], c=data['Q'], cmap='seismic',
                     marker='.')
    fig.colorbar(C1, ax=ax1)
    # Plot SI(Pr) v tlag
    C2 = ax2.scatter(data['FAST'], data['SI(Pr)'], c=data['Q'], cmap='seismic',
                     marker='.')
    fig.colorbar(C2, ax=ax2)
    ax2.set_xlabel('Fast direction')
    ax1.set_ylabel('SI(Pa)')
    ax2.set_ylabel('SI(Pr)')
    plt.show()

def plot_SI_diff_v_dt(data):
    '''Plots |SI(Pa) - SI(Pr)| v dt'''
    fig, ax = plt.subplots(1, 1, figsize=(7,7))
    data = data[data['Q']>0.5]
    dSI = abs(data['SI(Pa)'] - data['SI(Pr)'])
    
    C = ax.scatter(data['TLAG'], dSI, c=data['Q'], marker='.')
    fig.colorbar(C, label='Q')
    ax.set_ylabel('|SI(Pa) - SI(Pr)|')
    ax.set_xlabel('Delay time (s)')
    plt.show()
    
def scatter_bowtie(data_in):
    '''Plot a bowtie plot coloured by dt '''
    data = data_in[data_in['Q'] > 0.5]
    fig, ax = plt.subplots(1,1, figsize=(7,7))
    
    C = ax.scatter(data['SI(Pa)'], data['SI(Pr)'], c=data['TLAG'], marker='.')
    ax.plot([-5,0,5],[-5,0,5],'k--') # add y=x
    fig.colorbar(C, label='Lag Time')
    ax.set_xlabel('SI by approximation')
    ax.set_ylabel('SI by projection')
    plt.show()

if __name__ == '__main__':
    # Make sure we preserve eading 0s in time stamps by reading as strings
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}   
    sks = pd.read_csv(f'{RES_PATH}/Combined_SKS_sheba_results.sdb',
                      delim_whitespace=True,
                      converters=date_time_convert)
    skks = pd.read_csv(f'{RES_PATH}/Combined_SKKS_sheba_results.sdb',
                      delim_whitespace=True,
                      converters=date_time_convert)
    scs =  pd.read_csv(f'{RES_PATH}/Jacks_ScS_sheba_results.sdb',
                      delim_whitespace=True,
                      converters=date_time_convert)
    synth = read_and_cat_synth()
        
    SI_bowties(sks, skks, scs, synth)
    

    