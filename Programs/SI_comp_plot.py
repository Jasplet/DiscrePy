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

RES_PATH = '/Users/ja17375/DiscrePy/Sheba/Results/'

   
def read_and_cat_synth():
    
    df = pd.DataFrame()
    for spol in ['30', '45', '60']:
        for noise in ['01', '025', '005']:
             sp = pd.read_csv(f'{RES_PATH}/SYNTH/Noise{noise}/SP{spol}/SP{spol}_noise_{noise}_SYNTH_sheba_results.sdb',
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
    plt.savefig('/Users/ja17375/Documents/Thesis-enclosing/Thesis/chapters/chapter02/Figs/SI_bowties.png', dpi=600)


if __name__ == '__main__':
    # Make sure we preserve eading 0s in time stamps by reading as strings
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}   
    sks = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/Combined/Filt_05Hz/Combined_SKS_sheba_results.sdb',
                      delim_whitespace=True,
                      converters=date_time_convert)
    skks = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/Combined/Filt_05Hz/Combined_SKKS_sheba_results.sdb',
                      delim_whitespace=True,
                      converters=date_time_convert)
    scs =  pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/ScS/Jacks_ScS_sheba_results.sdb',
                      delim_whitespace=True,
                      converters=date_time_convert)
    synth = read_and_cat_synth()
        
    SI_bowties(sks, skks, scs, synth)


    