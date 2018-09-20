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
import misc

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

    def plot_fast_v_lag(self,ax,data,c='b',t=r'Plot of $\delta t$ versus $\phi$'):
        ''' Plots Fast direction versus lag time for a given set of events'''

        ax.errorbar(x=data.TLAG_SKS,y=data.FAST_SKS,yerr=data.DFAST_SKS,xerr=data.DTLAG_SKS,color='r',fmt='o',label='SKS')
        ax.errorbar(x=data.TLAG_SKKS,y=data.FAST_SKKS,yerr=data.DFAST_SKKS,xerr=data.DTLAG_SKKS,color='b',fmt='o',label='SKKS')
        for i in range(0,len(data.TLAG_SKS)):
            # print(i)
            ax.plot(x=[data.TLAG_SKS.values[i],data.TLAG_SKKS.values[i]],y=[data.FAST_SKS.values[i],data.FAST_SKKS.values[i]],color='black',ls='solid')
        ax.set_ylim([-90,90])
        ax.set_xlim([0,4])
        ax.set_xlabel(r'Lag Time $\delta t$ [s]')
        ax.set_ylabel(r'Fast Direction $\Phi$ [$\degree$]')
        ax.set_title(t)

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

        # ax1.legend(handles,loc='center left', bbox_to_anchor=(1.0, 0.5))
        ax1.set_title(r'$\phi$ v backazimuth') # $\delta$t values for the {:03d} SK(K)S pairs in bin {:04d}'.format(len(self.bin),self.bn))
        ax2.set_title(r'$\delta$t v backazimuth')

    def plot_baz_l2_dSI(self,ax1,ax2):
        ''' Make plot of Lambda 2 and d_SI v BAZ (of midpoints) for SKS and SKKS'''
        # fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize = (6,10))


        #plot_nulls
        ax1.plot(self.bin.Mid_BAZ,self.bin.LAM2,'.k',label='Lambda 2')
        ax2.plot(self.bin.Mid_BAZ,self.bin.D_SI,'.k',label='D_SI')


        lim = [np.round(np.min(self.bin.Mid_BAZ) - 5),np.round(np.max(self.bin.Mid_BAZ) + 5)]
        ax1.set_xlim(lim)
        ax1.set_ylim([0,1.0])
        ax1.set_ylabel(r'$\lambda _2$ values',fontsize=16)
        #
        #
        ax2.set_xlim(lim)
        ax2.set_ylim([0.,4.])
        ax2.set_ylabel(r'$\Delta _SI$ values',fontsize=16)
        ax2.set_xlabel(r'Backazimuth $[\degree]$',fontsize=16)

        # ax1.legend(handles,loc='center left', bbox_to_anchor=(1.0, 0.5))
        ax1.set_title(r'$\lambda_2$ v backazimuth',fontsize=18) # $\delta$t values for the {:03d} SK(K)S pairs in bin {:04d}'.format(len(self.bin),self.bn))
        ax2.set_title(r'$\Delta SI$ v backazimuth',fontsize=18)
        ax1.tick_params('both',labelsize=14)
        ax2.tick_params('both',labelsize=14)

    def plot_dSI(self,ax,SI,c='darkorange',t=r'Histogram of $\Delta$SI '):#save=False):
        '''plot a histogram dSI for the bin'''
        ax2 = ax.twinx() # Make a secondary y Axis
        # fig,ax = plt.subplots(1,1,figsize=(5,5))
        bins = np.arange(start=0,stop=3.2,step=0.2)
        h = ax.hist(SI,bins=bins,color=c,density=1)
        p1 = ax.plot([0.4,0.4],[0,20],color='black',ls='dashed',label='Dengs Threshold')
        p2 = ax.plot([self.avg_dSI(SI),self.avg_dSI(SI)],[0,max(h[0])],color='black',ls='solid',label=r'Median $\Delta SI$ = {:4.3f}'.format(self.avg_dSI(SI)))
        ylim = [0,max(h[0])]
        txt_y = max(h[0])/2.
        # Plot expected (truncated) normal distribution(s)
        x = np.linspace(0,3,1000)
        mu = np.mean(self.bin.D_SI)
        sigma = np.std(self.bin.D_SI)
        # print(mu,sigma)
        # d = [misc.std_norm(i - mu /sigma) for i in x]
        dt = [misc.trunc_norm(u=mu,s=sigma,x=i) for i in x]

        # ax.plot(x,d,linestyle='dashed')
        p3 = ax2.plot(x,dt,'b--',label = 'Expected Distribution')
        plots = p1+p2+p3
        labels = [l.get_label() for l in plots]
        ax2.set_ylim([0,0.5])
        ax.set_ylim(ylim)
        ax.set_xlim([0,max(h[1])])
        ax.set_ylabel('Density',fontsize=16)
        ax.set_xlabel(r'$\Delta$SI',fontsize=16)
        ax.set_title(t,fontsize=18)
        # ax.text(2,txt_y,r'Median $\Delta$SI = {:4.3f}'.format(self.avg_dSI(SI)))
        ax.legend(plots,labels,loc='best',fontsize=14)
        ax.tick_params('both',labelsize=14)

    def plot_lam2(self,ax,l2,c='darkorange',t=r'Histogram of $\lambda _2$ values'): #,save=False):
        '''plot a histogram of LAM2 values for the bin'''
        # fig,ax = plt.subplots(1,1,figsize=(6,6))
        ax2 = ax.twinx() # Make a secondaty y axis
        bins = np.arange(start=0,stop=1.1,step=0.05)
        h = ax.hist(l2,bins=bins,color=c,density=1)
        p1 = ax.plot([self.avg_lam2(l2),self.avg_lam2(l2)],[0,max(h[0])],color='black',ls='solid',label=r'Median $\lambda_2$ = {:4.3f}'.format(self.avg_lam2(l2)))
        txt_y = max(h[0])/2.
        ylim = [0,max(h[0])]
        # Plot expected (truncated) normal distributions
        x = np.linspace(0,3,1000)
        mu = np.mean(self.bin.LAM2)
        sigma = np.std(self.bin.LAM2)
        # print(mu,sigma)
        # d = [misc.std_norm(i - mu /sigma) for i in x]
        dt = [misc.trunc_norm(u=mu,s=sigma,x=i) for i in x]

        # ax.plot(x,d,linestyle='dashed')
        p2 = ax2.plot(x,dt,'--',color='blue',label='Expected Distribution')
        plots = p1+p2
        labels = [l.get_label() for l in plots]
        ax.set_ylim([0,0.5])
        ax.set_ylim(ylim)
        ax.set_xlim([0, max(h[1])])
        ax.set_ylabel('Density',fontsize=18)
        ax.set_xlabel(r'$\lambda _2$ value',fontsize=18)
        ax.set_title(t,fontsize=20)
        ax.legend(plots,labels,loc='best',fontsize=14)
        ax.tick_params('both',labelsize=14)
        # ax.text(0.7,txt_y,r'Median $\lambda_2$ = {:4.3f}'.format(self.avg_lam2(l2)))


    def avg_lam2(self,l2=None):
        '''return average lambda 2 value assuming a guassian distribution
            - lam2 : list/np array containing lambda 2 values

        '''
        # Currently just using numpy.median()
        if l2 is None:
            avg = np.median(self.bin.LAM2)
        else:
            avg = np.median(l2)
        return avg

    def avg_dSI(self,SI=None):
        '''return average lambda 2 value assuming a gaussian distribution in the bin
            - dSI  : list/np array containing Splitting intensity values
        '''
        if SI is None:
            avg = np.median(self.bin.D_SI)
        else:
            avg = np.median(SI)
        return avg

    def avg_splitting(self,fast,lag):
        '''return average FAST,LAG for the bin. Assuming a gaussian distribution
            - fast : list/np array of fast directions
            - lag : list/np array of lag times
        '''
        avg_f = np.median(self.bin.FAST_SKS)
        avg_l = np.median(self.bin.LAG_SKS)


    def plot(self,save=False):
        '''Make combined figure of BAZplots and Lam2/dSI histograms'''
        fig = plt.figure(figsize=(18,12))
        gs = gridspec.GridSpec(2,2)
        # ax = plt.subplot(gs[:,0:2])
        ax1 = plt.subplot(gs[0,0])
        ax2=plt.subplot(gs[1,0])
        ax3=plt.subplot(gs[0,1])
        ax4=plt.subplot(gs[1,1])
        # ax5=plt.subplot(gs[0,2])
        # ax6=plt.subplot(gs[1,2])
        self.plot_baz_l2_dSI(ax1,ax2)
        # self.plot_fast_v_lag(ax1,self.bin)

        self.plot_lam2(ax3,self.bin.LAM2)
        self.plot_dSI(ax4,self.bin.D_SI)

        # self.plot_lam2(ax5,self.bin[(self.bin.Q_SKS >= 0.5) | (self.bin.Q_SKKS >= 0.5)].LAM2,c='darkorange',t=r'$\lambda _2$ for Split SKS or SKKS')
        # self.plot_dSI(ax6,self.bin[(self.bin.Q_SKS >= 0.5) | (self.bin.Q_SKKS >= 0.5)].D_SI,c='darkorange',t=r'$\Delta$SI for Split SKS or SKKS')
        fig.suptitle(r'Analysis plots for trigonal bin no. {:04d} which contains {:03d} SK(K)S pairs'.format(self.bn,len(self.bin)),fontsize=20)
        # Either save the figure to the output directory or display it now
        if save is True:
            plt.savefig('{}/Analysis_plots_bin_{:04d}_{:03d}_pairs.eps'.format(self.fig_path,self.bn,len(self.bin)),format='eps',transparent=True)
            plt.close('all')
        elif save is False:
            plt.show()

def run(bins_file,plot,lim=10):
    ''' Function to run bin_analysis when imported in ipython environment'''

    bf = pd.read_csv(bins_file,delim_whitespace=True,converters={'TIME': lambda x: str(x),'DATE': lambda x : str(x)})

    counts = bf.bin_no.value_counts().copy() # make a dataeframe with the bin number and the count of how many times it occurs
    print('There are {} bins with {} pairs'.format(len(counts),len(bf)))
    print('Highest count is {} in bin {}'.format(counts[counts.idxmax()],counts.idxmax()))
    l2 = [ ] # Initialise lists to hold the average lambda 2 values for each bin
    dSI = [ ] # Initialise list to hold the average delta SI value for each bin
    no,cts,lat,long,V1_lat,V1_long,V2_lat,V2_long,V3_lat,V3_long = [ ],[ ],[ ],[ ],[ ],[ ],[ ],[ ],[ ],[ ]
    for i in counts.index:
        if counts[i] >= lim:
            # print(i)
            B = Bin(bf,bin_no=i)
            if plot is True:
                B.plot(save=True)

            no.append(i)
            cts.append(counts[i])
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
        else :
            print('Bin {} has {} counts, which is less than lim of {}'.format(i,counts[i],lim))

    # Make a dictioary of the series that we want to combine to make the datatframe
    dict = {'Bin_no' : no , 'Count' : cts , 'avg_lam2' : l2 ,
            'avg_dSI' : dSI , 'bin_lat' : lat, 'bin_long' : long , 'V1_lat' : V1_lat ,
            'V1_long' : V1_long , 'V2_lat' : V2_lat , 'V3_lat' : V3_lat , 'V3_long' : V3_long}

    # Make a dataframe for each bin

    df = pd.DataFrame(dict)
    return(df.round(4))


if __name__ == '__main__':
    print('Hello I am bin_analysis.py! You are running me form the Command line!')
    bins_file = sys.argv[1]
    if len(sys.argv) == 2:
        plot = False # Bool - switch for if you want to produce plots for the bins
    else:
        plot = sys.argv[2]

    run(bins_file,plot) # Call the run function
