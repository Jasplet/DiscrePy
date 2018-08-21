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
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


class Bin:
    '''A class to hold a signle bin and all in containing events. Requires pairs file (as dataframe) and desired Bin_number as inputs'''
    def __init__(self, df, bin_no,path='/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Combined/Bin_figs'):
        self.bn = bin_no # The bin number for the trigonal bin in questions
        self.bin = df[df.bin_no == bin_no].copy()
        self.fig_path = path # This path to the directory where the figures will be saved
        self.lat = self.bin.bin_lat.values[0]
        self.long = self.bin.bin_long.values[0]
        self.V1_lat = self.bin.V1_lat.values[0]
        self.V1_long = self.bin.V1_long.values[0]
        self.V2_lat = self.bin.V2_lat.values[0]
        self.V2_long = self.bin.V2_long.values[0]
        self.V3_lat = self.bin.V3_lat.values[0]
        self.V3_long = self.bin.V3_long.values[0]
        #parse nulls


    def plot(self,save=False):
        '''Make combined figure of BAZplots and Lam2/dSI histograms'''
        fig = plt.figure(figsize=(12,12))
        gs = gridspec.GridSpec(2,2)
        ax1 = plt.subplot(gs[0,0])
        ax2=plt.subplot(gs[1,0])
        ax3=plt.subplot(gs[0,1])
        ax4=plt.subplot(gs[1,1])
        self.plot_baz(ax1,ax2)
        self.plot_lam2(ax3)
        self.plot_dSI(ax4)

        fig.suptitle(r'Analysis plots for trigonal bin no. {:04d} which contains {:03d} SK(K)S pairs'.format(self.bn,len(self.bin)))
        # Either save the figure to the output directory or display it now
        if save is True:
            plt.savefig('{}/Analysis_plots_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
            plt.close('all')
        elif save is False:
            plt.show()


    def plot_baz(self,ax1,ax2):
        ''' Make plot of Fast and Lag v BAZ for SKS and SKKS'''
        # fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize = (6,10))

        #parse SKS nulls
        sks_null = self.bin[self.bin.Q_SKS < -0.7]
        skks_null = self.bin[self.bin.Q_SKKS < -0.7]
        #parse the clear split results
        sks_split = self.bin[self.bin.Q_SKS > 0.7]
        skks_split = self.bin[self.bin.Q_SKKS > 0.7]
        sks_meh = self.bin.drop(sks_null.index.append(sks_split.index)) # cuts out the clear splits and the clear nulls
        skks_meh = self.bin.drop(skks_null.index.append(skks_split.index))
        #plot the OK results
        # ax1.errorbar(x=sks_meh.BAZ,y=sks_meh.FAST_SKS,yerr=sks_meh.DFAST_SKS,color='blue',fmt='.',label='sks ok')
        # ax1.errorbar(x=skks_meh.BAZ,y=skks_meh.FAST_SKKS,yerr=skks_meh.DFAST_SKKS,color='blue',fmt='x',label='skks ok')
        #plot_nulls
        ax1.errorbar(x=sks_null.BAZ,y=sks_null.FAST_SKS,yerr=sks_null.DFAST_SKS,color='black',fmt='.',label='sks nulls')
        ax1.errorbar(x=skks_null.BAZ,y=skks_null.FAST_SKKS,yerr=skks_null.DFAST_SKKS,color='black',fmt='x',label='skks nulls')
        #plot_split
        ax1.errorbar(x=sks_split.BAZ,y=sks_split.FAST_SKS,yerr=sks_split.DFAST_SKS,color='darkorange',fmt='.',label='sks split')
        ax1.errorbar(x=skks_split.BAZ,y=skks_split.FAST_SKKS,yerr=skks_split.DFAST_SKKS,color='darkorange',fmt='x',label='skks split')

        lim = [np.round(np.min(self.bin.BAZ) - 5),np.round(np.max(self.bin.BAZ) + 5)]
        ax1.set_xlim(lim)
        ax1.set_ylim([-90,90])

        #plot the OK results
        # ax2.errorbar(x=sks_meh.BAZ,y=sks_meh.FAST_SKKS,yerr=sks_meh.DTLAG_SKS,color='blue',fmt='.',label='sks ok')
        # ax2.errorbar(x=skks_meh.BAZ,y=skks_meh.FAST_SKKS,yerr=skks_meh.DTLAG_SKKS,color='blue',fmt='x',label='skks ok')
        #plot_nulls
        ax2.errorbar(x=sks_null.BAZ,y=sks_null.TLAG_SKS,yerr=sks_null.DTLAG_SKS,color='black',fmt='.',label='sks nulls')
        ax2.errorbar(x=skks_null.BAZ,y=skks_null.TLAG_SKKS,yerr=skks_null.DTLAG_SKKS,color='black',fmt='x',label='skks nulls')
        #plot_split
        ax2.errorbar(x=sks_split.BAZ,y=sks_split.TLAG_SKS,yerr=sks_split.DTLAG_SKS,color='darkorange',fmt='.',label='sks split')
        ax2.errorbar(x=skks_split.BAZ,y=skks_split.TLAG_SKKS,yerr=skks_split.DTLAG_SKKS,color='darkorange',fmt='x',label='skks split')

        ax2.set_xlim(lim)
        ax2.set_ylim([0.,4.])

        ax1.legend(loc=0)
        ax1.set_title(r'$\phi$ v backazimuth') # $\delta$t values for the {:03d} SK(K)S pairs in bin {:04d}'.format(len(self.bin),self.bn))
        ax2.set_title(r'$\delta$t v backazimuth')

        # Either save the figure to the output directory or display it now
        # if save is True:
        #     plt.savefig('{}/BAZ_plot_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
        #     plt.close(fig)
        # elif save is False:
        #     plt.show()

    def plot_dSI(self,ax):#save=False):
        '''plot a histogram dSI for the bin'''

        # fig,ax = plt.subplots(1,1,figsize=(5,5))
        bins = np.arange(start=0,stop=3.2,step=0.2)
        h = ax.hist(self.bin.D_SI,bins=bins)
        ax.plot([0.4,0.4],[0,20],color='black',ls='dashed')
        ylim = [0,max(h[0])]
        ax.set_ylim(ylim)
        ax.set_ylabel('Frequency')
        ax.set_xlabel(r'$\Delta$SI')
        ax.set_title(r'Histogram of $\Delta$SI ')


        # Either save the figure to the output directory or display it now
        # if save is True:
        #     plt.savefig('{}/dSI_histogram_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
        #     plt.close(fig)
        # elif save is False:
        #     plt.show()

    def plot_lam2(self,ax): #,save=False):
        '''plot a histogram of LAM2 values for the bin'''
        # fig,ax = plt.subplots(1,1,figsize=(6,6))
        bins = np.arange(start=0,stop=1.1,step=0.05)
        h = ax.hist(self.bin.LAM2,bins=bins)
        ylim = [0,max(h[0])]
        ax.set_ylim(ylim)
        ax.set_ylabel('Frequency')
        ax.set_xlabel(r'$\lambda _2$ value')
        ax.set_title(r'Histogram of $\lambda _2$ values')



        # Either save the figure to the output directory or display it now
        # if save is True:
        #     plt.savefig('{}/lam2_histogram_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
        #     plt.close(fig)
        # elif save is False:
        #     plt.show()


    def avg_lam2(self):
        '''return average lambda 2 value assuming a guassian distribution
            - lam2 : list/np array containing lambda 2 values

        '''
        # Currently just using numpy.median()
        avg = np.median(self.bin.LAM2)
        return avg

    def avg_dSI(self):
        '''return average lambda 2 value assuming a gaussian distribution in the bin
            - dSI  : list/np array containing Splitting intensity values
        '''
        avg = np.median(self.bin.D_SI)
        return avg

    def avg_splitting(self,fast,lag):
        '''return average FAST,LAG for the bin. Assuming a gaussian distribution
            - fast : list/np array of fast directions
            - lag : list/np array of lag times
        '''
        avg_f = np.median(self.bin.FAST_SKS)
        avg_l = np.median(self.bin.LAG_SKS)

def run(bins_file):
    ''' Function to run bin_analysis when imported in ipython environment'''

    bf = pd.read_csv(bins_file,delim_whitespace=True,converters={'TIME': lambda x: str(x),'DATE': lambda x : str(x)})

    counts = bf.bin_no.value_counts().copy() # make a dataeframe with the bin number and the count of how many times it occurs
    print('There are {} bins with {} pairs'.format(len(counts),len(bf)))
    print('Highest count is {} in bin {}'.format(counts[counts.idxmax()],counts.idxmax()))
    l2 = [ ] # Initialise lists to hold the average lambda 2 values for each bin
    dSI = [ ] # Initialise list to hold the average delta SI value for each bin
    lat,long,V1_lat,V1_long,V2_lat,V2_long,V3_lat,V3_long = [ ],[ ],[ ],[ ],[ ],[ ],[ ],[ ]
    for i in counts.index:
        B = Bin(bf,bin_no=i)
        l2.append(B.avg_lam2())
        dSI.append(B.avg_dSI())
        lat.append(B.bin.bin_lat.values[0])
        long.append(B.bin.bin_long.values[0])
        V1_lat.append(B.bin.V1_lat.values[0])
        V1_long.append(B.bin.V1_long.values[0])
        V2_lat.append(B.bin.V2_lat.values[0])
        V2_long.append(B.bin.V2_long.values[0])
        V3_lat.append(B.bin.V3_lat.values[0])
        V3_long.append(B.bin.V3_long.values[0])

    # Make a dictioary of the series that we want to combine to make the datatframe
    dict = {'Bin_no' : counts.index , 'Count' : counts.values , 'avg_lam2' : l2 ,
            'avg_dSI' : dSI , 'bin_lat' : lat, 'bin_long' : long , 'V1_lat' : V1_lat ,
            'V1_long' : V1_long , 'V2_lat' : V2_lat , 'V3_lat' : V3_lat , 'V3_long' : V3_long}

    # Make a dataframe for each bin

    df = pd.DataFrame(dict)
    return(df)


if __name__ == '__main__':
    print('Hello I am bin_analysis.py! You are running me form the Command line!')
    bins_file = sys.argv[1]
    run(bins_file) # Call the run function
