#! /usr/bin/env python
######################
# bin_analysis.py
######################
# A program for the analysis of trigonal bins of SKS-SKKS pairs, created using James Wookey's program geogeom
# This program assumes that the binning has already been done (a function to run the binning may be aded later)
# and that the input file has matched each output bin number to the list of pairs

#####################################################################################################
# Imports
##############
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



class Bin:
    '''A class to hold a signle bin and all in containing events. Requires pairs file (as dataframe) and desired Bin_number as inputs'''
    def __init__(self, df, bin_no):
        self.bin = df[df.bin_no == bin_no].copy()

    def plot_baz(self):
        ''' Make plot of Fast and Lag v BAZ for SKS and SKKS'''
        fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize = (12,6))

        ax1.errorbar(x=self.bin.BAZ,y=self.bin.FAST_SKS,yerr=self.bin.DFAST_SKS,fmt='k.',label='sks')
        ax1.errorbar(x=self.bin.BAZ,y=self.bin.FAST_SKKS,yerr=self.bin.DFAST_SKKS,fmt='kx',label='skks')
        lim = [np.round(np.min(self.bin.BAZ) - 5),np.round(np.max(self.bin.BAZ) + 5)]
        ax1.set_xlim(lim)
        ax1.set_ylim([-90,90])

        ax2.errorbar(x=self.bin.BAZ,y=self.bin.TLAG_SKS,yerr=self.bin.DTLAG_SKS,fmt='k.')
        ax2.errorbar(x=self.bin.BAZ,y=self.bin.TLAG_SKKS,yerr=self.bin.DTLAG_SKKS,fmt='kx')
        ax2.set_ylim([0.,4.])

        ax1.legend(loc=0)

        plt.tight_layout()
        plt.show()

    def plot_dSI(self):
        '''plot a histogram dSI for the bin'''

    def plot_lam2(self):
        '''plot a histogram of LAM2 values for the bin'''

    def avg_lam2():
        '''return average lambda 2 value assuming a guassian distribution'''

    def avg_dSI():
        '''return average lambda 2 value assuming a gaussian distribution in the bin'''

    def avg_splitting():
        '''return average FAST,LAG for the bin. Assuming a gaussian distribution'''



if __name__ == '__main__':
  print('Hello I am bin_analysis.py! I am not currently deisgn ed to work from the command line. Please trying importing me into your script or using IPYTHON  :-) ')
