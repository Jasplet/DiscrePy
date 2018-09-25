#! /anaconda3/envs/splitwavepy/bin/python

# Welcome to discre.py. This script (or maybe module??) is for testing whether
# SKS SKKS pairs exhibit matching or discrepent splitting. This will be done
# by using the splitting measured for each phase to correct the splitting on
# the other.
### Imports
import numpy as np
import pandas as pd
import obspy as ob
import sys
import os
from os import path
#import splitwavepy as sw
import matplotlib.pyplot as plt
from stack import Stacker,plot_stack
from glob import glob
from datetime import datetime
from matplotlib import gridspec
# Maybe some others

###############################################################################
class Tester:

    def __init__(self,pr,path,mode='man'):
        '''
        Tester contains data and functions to Stack SKS and SKKS lambda2 surfaces in order to test for discrepancy

        pr - Full path to the .pairs (if stacking has not been done) file. If a .pairs file is passed in them discre.py wil automatically
             attempt to stack the LamR surfacesself.
             Alternatively a .stk file can be provided if you just want to do plotting of already existing stacks. Discre.py will
             assume that stacking has been done and will not attempt to do any stackingself.

        path - Path to the corresponding Runs directory where the .LamR files exist.
        '''
        self.path =path

        date_time_convert = {'TIME': lambda x: str(x).zfill(4),'DATE': lambda x : str(x)}
        # print(pr.split('.')[0])
        if pr.split('.')[-1] == 'pairs':
            # print('Hello')
            #File extention is .pairs, lam2 values dont exist so we need to run the stack
            self.pairs = pd.read_csv(pr,delim_whitespace=True,converters=date_time_convert)
            # Assuming that the main pairs filestem was read in, also read the matching and discepent pairs files
            self.lam2 = [ ]

            self.pair_stack(mode)
        elif pr.split('.')[-1] == 'stk':
            # File extention is .stk, lam2 values alreayd exist
            print('.stk files alreayd exist, reading')
            self.stk = pd.read_csv(pr,delim_whitespace=True,converters=date_time_convert)
            # pmatch = '{}_matches.stk'.format(pr.split('.')[0])
            # pdiff ='{}_diffs.stk'.format(pr.split('.')[0])
            # self.stk_matches = pd.read_csv(pmatch,delim_whitespace=True,converters=date_time_convert)
            # self.stk_diffs = pd.read_csv(pdiff,delim_whitespace=True,converters=date_time_convert)
            self.lam2 = self.stk.LAM2.values
        else:
            print('Unrecongised file type')

    def pair_stack(self,mode='man'):
        ''' Runs Stacker for all the desired pairs (a .pairs file)'''

        if mode =='man':
            ext ='lamR'
            print('man stacking')
            rd = input('Enter whick results directory you want the stacks saved under >')
            out = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}/Stacks'.format(rd)
            if os.path.isdir(out) is False:
                print('{} does not exist, creating'.format(out))
                os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}/Stacks'.format(rd))
        elif mode == 'sheba':
            ext ='lam2'
        # print('Running')

        for i,f in enumerate(self.pairs.DATE.values):
            # print('It {}, time is {} '.format(i,str(datetime.now())))
            # First get the right DATE,TIME and STATION
            date,time,stat = self.pairs.DATE[i], self.pairs.TIME[i], self.pairs.STAT[i]
            fstem = '{}_{}_{}'.format(stat,date,time)

            lam2_stem = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path,stat,fstem,ext))
            # print(lam2_stem)
            print('{}/{}/SKS/{}??_SKS.lam2'.format(self.path,stat,fstem))
            if len(lam2_stem) is not 0:
                # I.e if glob has managed to find the sks lam2 surface file
                sks_lam2 = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path,stat,fstem,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.{}'.format(self.path,stat,fstem,ext))[0]
                Stk = Stacker(sks_lam2,skks_lam2,fstem,out)
                if mode == 'man':
                    self.lam2.append(Stk.lam2)
                elif mode == 'sheba':
                    self.lam2.append(Stk.sol[-1])
            else:
                fstem2 = '{}_{}'.format(stat,date)
                print('fstem2')
                print('{}/{}/SKS/{}_*_SKS.{}'.format(self.path,stat,fstem2,ext))
                sks_lam2 = glob('{}/{}/SKS/{}_*_SKS.{}'.format(self.path,stat,fstem2,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}_*_SKKS.{}'.format(self.path,stat,fstem2,ext))[0]
                # Now for a sanity check
                if (len(sks_lam2) is not 0) or (len(skks_lam2) is not 0):
                    Stk = Stacker(sks_lam2,skks_lam2,fstem2,out)
                    if mode == 'man':
                        self.lam2.append(Stk.lam2)
                    elif mode == 'sheba':
                        self.lam2.append(Stk.sol[-1])
                else:
                    #print('lam2 surfaces cannot be found, skipping')
                    pass
    #        Now lets get the lambda 2 values
        print('Lam2 max: {} Lam2 min: {}'.format(max(self.lam2),min(self.lam2)))

        # self.lam2 = lam2

    def plot_lam2(self):
        print('Plotting')
        x = np.arange(0,len(self.lam2))
        # Sort lambda2 values (make sure they are in ascending order)
        fig,(ax0,ax1,ax2) = plt.subplots(1,3,figsize=(24,8),sharey=True)
        lam2 = self.lam2.copy()
        lam2.sort()
        m = self.stk_matches.LAM2.values.copy() # matching mabda 2 values
        m.sort()

        d =self.stk_diffs.LAM2.values.copy() #"discrepant" lambda 2 values
        d.sort()
        # print(x,y)
        ax0.plot(np.arange(0,len(m)),m,'k.')
        ax0.set_ylabel(r'$\lambda _2$ values')
        ax0.set_title('Matching lam2, {} pairs'.format(len(m)))

        ax1.plot(np.arange(0,len(d)),d,'k.')
        ax1.set_title('Discrepent lam2, {} pairs'.format(len(d)))
        t = self.path.split('/')[-1] # Gets the name of the Runs directory
        ax2.plot(x,lam2,'k.')
        # ax2.ylabel(r'$\lambda _2$ values')
        ax2.set_yticks(np.arange(0,1.4,step=0.2))
        ax2.set_ylim([0,1.4])
        ax2.set_title(r'$\lambda _2$ values for {} dataset, {} pairs'.format(t,len(lam2)))
        plt.tight_layout()
        plt.show()

    def hist_lam2(self):
        '''Plot a stacked histogram of lambda 2 values for "discrepent" and "matching" accoridng to classic +/- 2-sigma test'''
        bins = np.arange(0,1,0.1)
        lambs= [self.stk_matches.LAM2.values,self.stk_diffs.LAM2.values]

        fig,ax0 = plt.subplots(1,1,figsize=(8,8))

        #Stacked Hisogram
        ax0.hist(lambs,bins, histtype='bar', stacked=True,label=['Matching','Discrepent'])
        ax0.legend(loc=0)
        ax0.set_title(r'$\lambda _2$ values for matching and discrepent splitting')
        ax0.set_ylabel('Frequency')
        ax0.set_xlabel(r'$\lambda _2$ values')
        plt.show()

    def lam2_v_SI(self,figname):
        '''Plot a scatter plot of lambda2 values verus splitting INtensity difference'''

        fig,ax = plt.subplots(1,1,figsize=(8,8))

        ax.plot(self.stk_matches.LAM2,self.stk_matches.D_SI,color='blue',marker='.',ls='None',label="'Matching'")
        ax.plot(self.stk_diffs.LAM2,self.stk_diffs.D_SI,color='darkorange',marker='.',ls='None',label="'Discrepant'")
        ax.plot([0, 1.0],[0.4,0.4],ls='dashed',color='black')
        ax.set_xlabel(r'$\lambda _2$ values')
        ax.set_ylabel(r'$\Delta$ SI')
        ax.legend()
        ax.set_xlim([0,1])
        ax.set_ylim([0,4.0])
        plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/{}.eps'.format(figname),format='eps',dpi=1000)
        plt.show()

    def discrepancy_plot(self,nplots=2,surfs_to_plot=None,save=False,sigma=1,**kwargs):
        '''Top level plotting function for surfaces to look for discrepancy in Splitting
            nplots - the number of plots that you want (if the surfs_to_plot is not specified)
            surfs_to_plot - allows
            sigma - multiplier to error bounds of splitting results
        '''
        self.p_sorted = self.stk.sort_values(by='LAM2',ascending=True)
        self.p_sorted.reset_index(drop=True)

        # Find indecies of events we want to plot
        if surfs_to_plot is None:
            self.surfs= np.round(np.arange(0,len(self.p_sorted),round((len(self.p_sorted)/nplots))))
        else:
            self.surfs = list(set([i if i < len(self.p_sorted) else (len(self.p_sorted)-1) for i in surfs_to_plot])) # This makes sure that indicies are always within the range of available surfaces (stops errors for occuring)
            self.surfs = self.surfs.sort() # Sorts list in ascending order, has to be done speratly as sort acts of list and returns nothing
        # print(self.surfs)
        if save is True:
            dir = input('Enter Directory you want to save stacked surfaces to > ')

            if os.path.isdir('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}'.format(dir)) is False:
                make_d = input('Directory {} does not exist, do you want to create it? (y/n)'.format(dir))
                if make_d =='y':
                    print('Ok, creating directory {}'.format(dir))
                    os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}'.format(dir))
                else:
                    print('Exiting....')
                    sys.exit()

        for s in self.surfs:
            stat,date,time = self.p_sorted.STAT.values[s],self.p_sorted.DATE.values[s], self.p_sorted.TIME.values[s]
            # Note that dtlag and dfast are multiplied through by sigma HERE !!
            fast_sks,dfast_sks,lag_sks,dlag_sks = self.p_sorted.FAST_SKS.values[s],(sigma*self.p_sorted.DFAST_SKS.values[s]),self.p_sorted.TLAG_SKS.values[s],(sigma*self.p_sorted.DTLAG_SKS.values[s])
            fast_skks,dfast_skks,lag_skks,dlag_skks = self.p_sorted.FAST_SKKS.values[s],(sigma*self.p_sorted.DFAST_SKKS.values[s]),self.p_sorted.TLAG_SKKS.values[s],(sigma*self.p_sorted.DTLAG_SKKS.values[s])
            lam2 = self.p_sorted.LAM2.values[s]
            print('Stat {}, Evt Time {}-{} LAM2 = {}'.format(stat,date,time,lam2))
            l_path = '{}/{}/{}_{}_{}'.format(self.path,stat,stat,date,time) #Path to lambda 2 surfaces for SKS and SKKS
            self.lam2_surface(l_path)

            # fig = plt.figure(figsize=(12,12))
            fig, (ax0,ax1,ax2) = plt.subplots(1,3,figsize=(24,8),sharey=True)
            fig.patch.set_facecolor('None')
            plt.suptitle(r'Event {}_{}_{}. Stacked $\lambda _2$ value = {:4.3f}'.format(stat,date,time,self.p_sorted.LAM2.values[s]),fontsize=28)
            # gs = gridspec.GridSpec(3,2)
            # ax0 = plt.subplot(gs[0,0])
            ax0.set_title(r'SKS $\lambda _2$ surface',fontsize=24)
            C0 = ax0.contourf(self.T,self.F,(self.sks_lam2),[0.005,0.01,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5],cmap='inferno_r',extend='both')
            ax0.contour(C0,colors='k')
            # ax0.clabel(C0,C0.levels,inline=True,fmt ='%2.3f')
            #Plot SKS Solution
            ax0.plot(lag_sks,fast_sks,'b.',label='SKS Solution')
            ax0.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax0.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            ax0.set_ylabel(r'Fast,$\phi$, (deg)')
            ax0.set_xlabel(r'Lag ,$\delta$ t, (sec)')
            #Plot SKKS Solution
            ax0.plot(lag_skks,fast_skks,'r.',label='SKKS Solution')
            ax0.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'r-')
            ax0.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'r-')
            ax0.set_ylim([-90,90])
            ax0.set_xlim([0,4])
            ax0.set_yticks([-90,-60,-30,0,30,60,90])
            # ax0.contourf(self.sks_lam2,cmap='inferno_r')
            # ax1 = plt.subplot(gs[0,1])
            C1 = ax1.contourf(self.T,self.F,self.skks_lam2,[0,0.005,0.01,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5],cmap='inferno_r',extend='both')
            C3 = ax1.contour(C1,colors='k')
            # ax1.clabel(C1,C1.levels,inline=True,fmt ='%2.3f')
            # ax1.contourf(self.skks_lam2,cmap='magma')
            #Plot SKS Solution
            ax1.plot(lag_sks,fast_sks,'b.',label='SKS Solution')
            ax1.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax1.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            #Plot SKKS Solution
            ax1.plot(lag_skks,fast_skks,'r.',label='SKKS Solution')
            ax1.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'r-')
            ax1.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'r-')
            ax1.set_xlabel(r'Lag ,$\delta$ t, (sec)')
            ax1.set_ylim([-90,90])
            ax1.set_xlim([0,4])
            ax1.set_yticks([-90,-60,-30,0,30,60,90])
            ax1.set_title(r'SKKS $\lambda _2$ surface',fontsize=24)
            # ax2 = plt.subplot(gs[1:,:])

            self.show_stacks(ax2,l_path.split('/')[-1])
            # print('STK_FAST: {} +/- {}'.format(self.stk_fast,self.stk_dfast))
            # Modify stk_dlag and stk_dfast by sigma
            ##########################
            # self.stk_dlag = self.stk_dlag*sigma
            # self.stk_dfast = self.stk_dfast*sigma
            ############################
            ax2.set_ylim([-90,90])
            ax2.set_xlim([0,4])
            ax2.set_yticks([-90,-60,-30,0,30,60,90])
            #Plto SKS solution
            ax2.plot(lag_sks,fast_sks,'b.',label='SKS Solution')
            ax2.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax2.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            #Plot SKKS solution
            ax2.plot(lag_skks,fast_skks,'r.',label='SKKS Solution')
            ax2.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'r-')
            ax2.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'r-')
            # Plot Stacked solution on SKS
            ax0.plot(self.stk_lag,self.stk_fast,'g.',label='Stacked Solution')
            # ax0.plot([self.stk_lag- self.stk_dlag, self.stk_lag + self.stk_dlag],[ self.stk_fast, self.stk_fast],'g-')
            # ax0.plot([self.stk_lag, self.stk_lag],[ self.stk_fast - self.stk_dfast ,self.stk_fast +self.stk_dfast],'g-')
            # Plot Stacked Solution on SKKS surface
            ax1.plot(self.stk_lag,self.stk_fast,'g.',label='Stacked Solution')
            # ax1.plot([self.stk_lag-self.stk_dlag,self.stk_lag+self.stk_dlag],[self.stk_fast,self.stk_fast],'g-')
            # ax1.plot([self.stk_lag,self.stk_lag],[self.stk_fast-self.stk_dfast,self.stk_fast+self.stk_dfast],'g-')
            ## Add a legend (on ax0)
            ax0.legend(bbox_to_anchor=(0,1),loc='upper left')


            # cb = fig.colorbar(C1)
            # cb.add_lines(C3)

            ax2.set_title('Stacked SKS SKKS surface',fontsize=24)
            if save is True:
                # dir = input('Enter Directory you want to save stacked surfaces to > ')
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}/LAM2_{:4.4f}_STAT_{}.eps'.format(dir,lam2,stat),format='eps',dpi=800)
                plt.close()


        # Read Lam2 surfaces for SKS and SKKS
        # self.lam2_surface(stem)

            if save is False:
                plt.show()

    def write_lam2(self,out_stem):
        '''Adds lam2 values to pairs'''
        self.stk =self.pairs.copy() # Create copy of pairs dataframe
        l2df = {'LAM2' : self.lam2}
        ldf = pd.DataFrame(l2df)
        self.stk['LAM2'] = ldf # Add Lambda 2 vlaues as a  new column to the dataframe
        outfile = '{}.stk'.format(out_stem.split('.')[0]) # Change file extension from .pairs to .stk to denote that the stacked lambda 2 value have been added
        print(outfile)
        self.stk.to_csv(outfile,sep=' ',index=False)

    def return_lam2(self):
        ''' Returns lambda2 values as a DF'''

        l2df = {'LAM2' : self.lam2} # Dict to create Datatframe

        return pd.DataFrame(l2df)  


    def lam2_surface(self,fstem):
        ''' Function to read  SKS and SKKS .lam2 surface files from sheba '''
        # print(fstem)
        t_sks = '{}/SKS/{}??_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),fstem.split('/')[-1])
        t_skks = '{}/SKKS/{}??_SKKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),fstem.split('/')[-1])
        sks =glob(t_sks)
        skks = glob(t_skks)
        # print(t_sks)
        if len(sks) == 0:
            stem = '_'.join(fstem.split('/')[-1].split('_')[0:-1]) # aka cut off time part of file extension
            # print('{}/SKS/{}*_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),fstem.split('/')[-1]))
            sks = glob('{}/SKS/{}*_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),stem))
            skks = glob('{}/SKKS/{}*_SKKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),stem))

        self.sks_lam2 = np.loadtxt(sks[0])#,skiprows=4) # skip rows not needed for .lamR files
        self.skks_lam2 = np.loadtxt(skks[0])#,skiprows=4)

        nfast,nlag = self.sks_lam2.shape ;
        lag_max = 4.
        [self.T,self.F] = np.meshgrid(np.linspace(0,lag_max,num=nlag),np.arange(-90,91,1)) ;

    def show_stacks_sheba(self,ax,path,evt):
        '''Function to find and plot desired surface stacks based on the LAMDA2 value '''
        ### Plot Min Lamnda 2

        with open('{}/sheba_stack.sol'.format(path),'r') as reader:
            head = reader.readline()  #Reads headers
            S = reader.readline().split() # Reads solution
            fast,dfast = float(S[0]), float(S[1])
            lag,dlag = float(S[2]),float(S[3])
            nsurf = float(S[4])
            lag_step = float(S[5])
            lam2 = S[6]
            # print(lam2)
        # Read surface
        err = np.loadtxt('{}/sheba_stack.err'.format(path))
        nfast,nlag = err.shape ;
        lag_max = (nlag) * lag_step ;
        [T,F] = np.meshgrid(np.arange(0,lag_max,lag_step),np.arange(-90,91,1)) ;

        C = ax.contour(T,F,err,15,colors='k')
        ax.set_ylabel(r'Fast,$\phi$, (deg)')
        ax.set_xlabel(r'Lag ,$\delta$ t, (sec)')
        ax.plot([lag-dlag,lag+dlag],[fast,fast],'g-')
        ax.plot([lag,lag],[fast-dfast,fast+dfast],'g-')
        ax.plot(lag,fast,'g.')
        ax.clabel(C,C.levels,inline=True,fmt ='%4.3f')
        # ax.set_title(r'Event {}. $\lambda$ 2 value = {}'.format(evt,lam2))
        # Add fast, lag as attributes so we can plot them elsewhere
        self.stk_fast,self.stk_dfast
        self.stk_lag,self.stk_dlag

    def show_stacks(self,ax,evt,path='/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Nulls'):
        '''Function to find and plot desired surface stacks based on the LAMDA2 value '''
        ### Plot Min Lamnda 2
        print('{}/Stacks/{}??.lamSTK'.format(path,evt))
        file =glob('{}/Stacks/{}??.lamSTK'.format(path,evt))
        if len(file) == 0:
            file = glob('{}/Stacks/{}*.lamSTK'.format(path,'_'.join(evt.split('_')[0:-1])))

        stk = np.loadtxt(file[0])
        nfast,nlag = stk.shape ;
        lag_step = 0.025
        lag_max = (nlag) * lag_step;
        [T,F] = np.meshgrid(np.arange(0,lag_max,lag_step),np.arange(-90,91,1)) ;
        jf,jt  = np.unravel_index(stk.argmin(),stk.shape)
        print(stk.argmin())
        print(jf,jt)

        fast = np.arange(-90,91,1)[jf]
        lag = np.arange(0,lag_max,lag_step)[jt]
        C = ax.contourf(T,F,stk,[0.01,0.025,0.05,0.075,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.5],cmap='inferno_r',extend='both')
        C2 = ax.contour(C,colors='k')
        ax.set_ylabel(r'Fast,$\phi$, (deg)')
        ax.set_xlabel(r'Lag ,$\delta$ t, (sec)')
        # ax.plot([lag-dlag,lag+dlag],[fast,fast],'g-')
        # ax.plot([lag,lag],[fast-dfast,fast+dfast],'g-')
        ax.plot(lag,fast,'g.')
        # ax.clabel(C,C.levels,inline=True,fmt ='%4.3f')

        # self.cbar = plt.colorbar(C)
        # self.cbar.add_lines(C2)
        # ax.set_title(r'Event {}. $\lambda$ 2 value = {}'.format(evt,lam2))
        # Add fast, lag as attributes so we can plot them elsewhere
        sol = stk.min()

        print('Lam2 {}, fast {} lag {}'.format(sol,fast,lag))
        self.stk_fast = fast
        self.stk_lag = lag

if __name__ == '__main__':
    print('This is discre.py')
    # p = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Accepted_SKS_SKKS_all.pairs'
    # path = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split'
    p = sys.argv[1]#'/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Deng/Deng_events_SKS_SKKS.pairs'
    path = sys.argv[2]#'/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Deng_events'
    print('File {} Path {}'.format(p,path))
    t = Tester(p,path,mode='man')

    t.write_lam2(p) # Writes lam2 values to new .stk file

    # show_stacks(p2)
    #print(lam2)
    #plot_lam2(p.index.values,lam2)
