#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 14:56:16 2021

@author: ja17375

Reproduce SNR test figure for synthetics
"""
import pandas as pd
import matplotlib.pyplot as plt
DATA_PATH = '~/DiscrePy/Sheba/Results/SYNTH/SNR_test'

def plot_snr_v_splitting(synth):

    fig, axs = plt.subplots(2, 2, sharex=True, figsize=(7,6))
    
    #plot fast dir
    axs[0,0].plot(synth.SNR, synth.FAST, 'k.')
    #plot fast dir error
    axs[1,0].plot(synth.SNR, synth.DFAST, 'k.')
    axs[0,0].set_ylabel(r'Fast direction, $\phi$ ($\degree$)')
    axs[0,0].set_ylim([-90,90])
    axs[1,0].set_ylabel(r'Fast direction error, $\sigma_{\phi}$ ($\degree$)')
    axs[1,0].set_ylim([0,45])
    axs[1,0].set_xlabel('Signal/Noise ratio')
    axs[1,0].set_xlim([0,20])
    
    axs[0,1].plot(synth.SNR, synth.TLAG, 'k.')
    axs[0,1].set_ylim([0,4])
    axs[0,1].set_ylabel('Lag time, $\delta t$ (s)')
    axs[1,1].plot(synth.SNR, synth.DTLAG, 'k.')
    axs[1,1].set_ylim([0,1.5])
    axs[1,1].set_xlabel('Signal/Noise ratio')
    axs[1,1].set_ylabel('Lag time error, $\sigma_{\delta t}$ (s)')
    
    # Add SNR 5 cutoff
    axs[0,0].axvline(5, ymin=-90, ymax=90, linestyle='dashed', color='black')
    axs[1,0].axvline(5, ymin=0, ymax=45, linestyle='dashed', color='black')
    axs[0,1].axvline(5, ymin=0, ymax=4, linestyle='dashed', color='black')
    axs[1,1].axvline(5, ymin=0, ymax=1.5, linestyle='dashed', color='black')
    
    return fig
if __name__ == '__main__':
    
    S1= pd.read_csv(f'{DATA_PATH}/SP30_SNR_range_test_P1_SYNTH_sheba_results.sdb', delim_whitespace=True)
    S2= pd.read_csv(f'{DATA_PATH}/SP30_SNR_range_test_P2_SYNTH_sheba_results.sdb', delim_whitespace=True)
    S3= pd.read_csv(f'{DATA_PATH}/SP30_SNR_range_test_P3_SYNTH_sheba_results.sdb', delim_whitespace=True)
#   S1 (phi =-10, dt = 2.0s)
#   S2 (phi = 70, dt = 1.5s)
#   S3 (phi = 50, dt = 2.5s)
    fig = plot_snr_v_splitting(S1)
    fig.suptitle(r'Synthetics: $\phi = -10\degree$, $\delta t = 2.0$s')
    fig.tight_layout()
    fig.savefig('/Users/ja17375/Documents/Thesis-enclosing/Thesis/chapters/chapter02/Figs/Synth_SNR_S1.png', dpi=600)
    
    fig2 = plot_snr_v_splitting(S1)
    fig2.suptitle(r'Synthetics: $\phi = 70\degree$, $\delta t = 1.5$s')
    fig2.tight_layout()
    fig2.savefig('/Users/ja17375/Documents/Thesis-enclosing/Thesis/chapters/chapter02/Figs/Synth_SNR_S2.png', dpi=600)

    fig3 = plot_snr_v_splitting(S1)
    fig3.suptitle(r'Synthetics: $\phi = 70\degree$, $\delta t = 1.5$s')
    fig3.tight_layout()
    fig3.savefig('/Users/ja17375/Documents/Thesis-enclosing/Thesis/chapters/chapter02/Figs/Synth_SNR_S3.png', dpi=600)
