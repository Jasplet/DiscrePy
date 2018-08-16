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
    def __init__(self, df, bin_no,path='/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Combined/Bin_figs'):
        self.bn = bin_no # The bin number for the trigonal bin in questions
        self.bin = df[df.bin_no == bin_no].copy()
        self.fig_path = path # This path to the directory where the figures will be saved

    def plot_baz(self,save=False):
        ''' Make plot of Fast and Lag v BAZ for SKS and SKKS'''
        fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize = (6,10))

        ax1.errorbar(x=self.bin.BAZ,y=self.bin.FAST_SKS,yerr=self.bin.DFAST_SKS,fmt='k.',label='sks')
        ax1.errorbar(x=self.bin.BAZ,y=self.bin.FAST_SKKS,yerr=self.bin.DFAST_SKKS,fmt='kx',label='skks')
        lim = [np.round(np.min(self.bin.BAZ) - 5),np.round(np.max(self.bin.BAZ) + 5)]
        ax1.set_xlim(lim)
        ax1.set_ylim([-90,90])

        ax2.errorbar(x=self.bin.BAZ,y=self.bin.TLAG_SKS,yerr=self.bin.DTLAG_SKS,fmt='k.')
        ax2.errorbar(x=self.bin.BAZ,y=self.bin.TLAG_SKKS,yerr=self.bin.DTLAG_SKKS,fmt='kx')
        ax2.set_ylim([0.,4.])

        ax1.legend(loc=0)
        ax1.set_title(r'$\phi$ and $\delta$t values for the {:03d} SK(K)S pairs in bin {:04d}'.format(len(self.bin),self.bn))
        plt.tight_layout()
        # Either save the figure to the output directory or display it now
        if save is True:
            plt.savefig('{}/BAZ_plot_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
            plt.close(fig)
        elif save is False:
            plt.show()

    def plot_dSI(self,save=False):
        '''plot a histogram dSI for the bin'''

        fig,ax = plt.subplots(1,1,figsize=(5,5))
        bins = np.arange(start=0,stop=4.2,step=0.2)
        h = ax.hist(self.bin.D_SI,bins=bins)

        ax.set_ylabel('Frequency')
        ax.set_xlabel(r'$\Delta$SI')
        ax.set_title(r'Histogram of $\Delta$SI for the {:03d} SK(K)S pairs in bin {:04d}'.format(len(self.bin),self.bn))
        plt.tight_layout()

        # Either save the figure to the output directory or display it now
        if save is True:
            plt.savefig('{}/dSI_histogram_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
            plt.close(fig)
        elif save is False:
            plt.show()

    def plot_lam2(self,save=False):
        '''plot a histogram of LAM2 values for the bin'''
        fig,ax = plt.subplots(1,1,figsize=(6,6))
        bins = np.arange(start=0,stop=1.1,step=0.1)
        h = ax.hist(self.bin.LAM2,bins=bins)
        ax.set_ylabel('Frequency')
        ax.set_xlabel(r'$\lambda _2$ value')
        ax.set_title(r'Histogram of $\lambda _2$ values for the {:03d} SK(K)S pairs in bin {:04d}'.format(len(self.bin),self.bn))
        plt.tight_layout

        # Either save the figure to the output directory or display it now
        if save is True:
            plt.savefig('{}/lam2_histogram_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
            plt.close(fig)
        elif save is False:
            plt.show()


    def avg_lam2():
        '''return average lambda 2 value assuming a guassian distribution'''

    def avg_dSI():
        '''return average lambda 2 value assuming a gaussian distribution in the bin'''

    def avg_splitting():
        '''return average FAST,LAG for the bin. Assuming a gaussian distribution'''



if __name__ == '__main__':
  print('Hello I am bin_analysis.py! I am not currently deisgn ed to work from the command line. Please trying importing me into your script or using IPYTHON  :-) ')
