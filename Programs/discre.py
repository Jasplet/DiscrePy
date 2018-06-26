#! /anaconda3/envs/splitwavepy/bin/python

# Welcome to discre.py. This script (or maybe module??) is for testing whether
# SKS SKKS pairs exhibit matching or discrepent splitting. This will be done
# by using the splitting measured for each phase to correct the splitting on
# the other.
### Imports
import numpy as np
import pandas as pd
import obspy as ob
#import splitwavepy as sw
import matplotlib.pyplot as plt
from stack import Stacker,plot_stack
from glob import glob
from datetime import datetime
from matplotlib import gridspec
# Maybe some others

###############################################################################
class Tester:

    def __init__(self,pr,path,stack_done=False):


        date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        if stack_done is False:
            self.pairs = pd.read_csv(pr,delim_whitespace=True,converters=date_time_convert)
            self.lam2 = [ ]
        elif stack_done is True:
            self.p2 = pd.read_csv(pr,delim_whitespace=True,converters=date_time_convert)
            self.lam2 = self.p2.LAM2.values

        self.path =path



    def pair_stack(self):
        ''' Runs Stacker for all the desired pairs (a .pairs file)'''

        print('Running')
        for i,f in enumerate(self.pairs.DATE.values):
            #rint('It {}, time is {} '.format(i,str(datetime.now())))
            # First get the right DATE,TIME and STATION
            date,time,stat = self.pairs.DATE[i], self.pairs.TIME[i], self.pairs.STAT[i]
            fstem = '{}_{}_{}'.format(stat,date,time)

            lam2_stem = glob('{}/{}/SKS/{}??_SKS.lam2'.format(self.path,stat,fstem))
            print(lam2_stem)
            # print('{}/{}/SKS/{}??_SKS.lam2'.format(self.path,stat,fstem))
            if len(lam2_stem) is not 0:
                # I.e if glob has managed to find the sks lam2 surface file
                sks_lam2 = glob('{}/{}/SKS/{}??_SKS.lam2'.format(self.path,stat,fstem))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.lam2'.format(self.path,stat,fstem))[0]
                Stk = Stacker(sks_lam2,skks_lam2,fstem)
                self.lam2.append(Stk.sol[-1])
            else:
                fstem2 = '{}_{}'.format(stat,date)
                sks_lam2 = glob('{}/{}_*_SKS.lam2'.format(self.path,fstem2))[0]
                skks_lam2 = glob('{}/{}_*_SKS.lam2'.format(self.path,fstem2))[0]
                # Now for a sanity check
                if (len(sks_lam2) is not 0) or (len(skks_lam2) is not 0):
                    Stk = Stacker(sks_lam2,skks_lam2)
                    self.lam2.append(Stk.sol[-1])
                else:
                    #print('lam2 surfaces cannot be found, skipping')
                    pass
    #        Now lets get the lambda 2 values
        print('Lam2 max: {} Lam2 min: {}'.format(max(self.lam2),min(self.lam2)))

        # self.lam2 = lam2

    def plot_lam2(self):
        print('Plotting')

        plt.plot(np.arange(0,len(self.lam2)),self.lam2.sort(),'k.')
        plt.ylabel('lambda 2 values')
        plt.yticks(np.arange(0,2,step=0.2))
        plt.ylim([0,2])
        plt.show()

    def discrepancy_plot(self,nplots=2,surfs_to_plot=None,save=False,**kwargs):
        '''Top level plotting function for surfaces to look for discrepancy in Splitting
            nplots - the number of plots that you want (if the surfs_to_plot is not specified)
            surfs_to_plot - allows
        '''
        self.p_sorted = self.p2.sort_values(by='LAM2',ascending=True)
        # print(self.p_sorted)
        # Find indecies of events we want to plot
        if surfs_to_plot is None:
            self.surfs= np.round(np.arange(0,len(self.p_sorted),round((len(self.p_sorted)/nplots))))
        else:
            self.surfs = list(set([i if i < len(self.p_sorted) else (len(self.p_sorted)-1) for i in surfs_to_plot])) # This makes sure that indicies are always within the range of available surfaces (stops errors for occuring)
            self.surfs.sort() # Sorts list in ascending order, has to be done speratly as sort acts of list and returns nothing
        print(self.surfs)
        if save is True:
            dir = input('Enter Directory you want to save stacked surfaces to > ')

        for s in self.surfs:
            # print(s)
            stat,date,time = self.p_sorted.STAT.values[s],self.p_sorted.DATE.values[s], self.p_sorted.TIME.values[s]
            fast_sks,dfast_sks,lag_sks,dlag_sks = self.p_sorted.FAST_SKS.values[s],self.p_sorted.DFAST_SKS.values[s],self.p_sorted.TLAG_SKS.values[s],self.p_sorted.DTLAG_SKS.values[s]
            fast_skks,dfast_skks,lag_skks,dlag_skks = self.p_sorted.FAST_SKKS.values[s],self.p_sorted.DFAST_SKKS.values[s],self.p_sorted.TLAG_SKKS.values[s],self.p_sorted.DTLAG_SKKS.values[s]
            lam2 = self.p_sorted.LAM2.values[s]
            print('Stat {}, Evt Time {}-{} LAM2 = {}'.format(stat,date,time,lam2))
            l_path = '{}/{}/{}_{}_{}'.format(self.path,stat,stat,date,time)
            print(l_path)
            self.lam2_surface(l_path)

            # fig = plt.figure(figsize=(12,12))
            fig, (ax0,ax1,ax2) = plt.subplots(1,3,figsize=(24,8),sharey=True)
            plt.suptitle(r'Event {}_{}_{}. Stacked $\lambda _2$ value = {}'.format(stat,date,time,self.p_sorted.LAM2.values[s]),fontsize=16)
            # gs = gridspec.GridSpec(3,2)
            # ax0 = plt.subplot(gs[0,0])
            ax0.set_title(r'SKS $\lambda _2$ surfaces')
            C0 = ax0.contour(self.T,self.F,self.sks_lam2,[0.5,1,2,3,4,5,10,15,20,25,30,40,50],colors='k')
            ax0.clabel(C0,C0.levels,inline=True,fmt ='%2.0f')
            ax0.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax0.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            ax0.set_ylim([-90,90])
            ax0.set_xlim([0,4])
            ax0.set_yticks([-90,-60,-30,0,30,60,90])
            # ax0.contourf(self.sks_lam2,cmap='magma_r')
            # ax1 = plt.subplot(gs[0,1])
            C1 = ax1.contour(self.T,self.F,self.skks_lam2,[0.5,1,2,3,4,5,10,15,20,25,30,40,50],colors='k')
            ax1.clabel(C1,C1.levels,inline=True,fmt ='%2.0f')
            # ax1.contourf(self.skks_lam2,cmap='magma')
            ax1.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'b-')
            ax1.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'b-')
            ax1.set_ylim([-90,90])
            ax1.set_xlim([0,4])
            ax1.set_yticks([-90,-60,-30,0,30,60,90])
            ax1.set_title(r'SKKS $\lambda _2$ surfaces')
            # ax2 = plt.subplot(gs[1:,:])

            self.show_stacks(ax2,'{}/'.format(l_path),l_path.split('/')[-1])
            ax2.set_ylim([-90,90])
            ax2.set_xlim([0,4])
            ax2.set_yticks([-90,-60,-30,0,30,60,90])
            plt.title('Stacked SKS SKKS surface')
            if save is True:
                # dir = input('Enter Directory you want to save stacked surfaces to > ')
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}/LAM2_{}_STAT_{}.png'.format(dir,lam2,stat))
                plt.close()


        # Read Lam2 surfaces for SKS and SKKS
        # self.lam2_surface(stem)

            if save is False:
                plt.show()

    def write_lam2(self):
        '''Adds lam2 values to pairs'''
        l2df = {'LAM2' : self.lam2}
        ldf = pd.DataFrame(l2df)
        self.pairs['LAM2'] = ldf

        self.pairs.to_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Accepted_SKS_SKKS_all_w_lam2.pairs',sep=' ')

        self.p2 = self.pairs

    def lam2_surface(self,fstem):
        ''' Function to read  SKS and SKKS .lam2 surface files from sheba '''
        sks =glob('{}/{}??_SKS.lam2'.format(fstem,fstem.split('/')[-1]))
        skks = glob('{}/{}??_SKKS.lam2'.format(fstem,fstem.split('/')[-1]))
        # print(sks)
        self.sks_lam2 = np.loadtxt(sks[0],skiprows=4)
        self.skks_lam2 = np.loadtxt(skks[0],skiprows=4)

        nfast,nlag = self.sks_lam2.shape ;
        lag_max = 4.
        [self.T,self.F] = np.meshgrid(np.linspace(0,lag_max,num=nlag),np.arange(-90,91,1)) ;

    def show_stacks(self,ax,path,evt):
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

        C = ax.contour(T,F,err,[1,2,3,4,5,10,15,20,50,100],colors='k')
        ax.set_ylabel(r'Fast,$\phi$, (deg)')
        ax.set_xlabel(r'Lag ,$\delta$ t, (sec)')
        ax.plot([lag-dlag,lag+dlag],[fast,fast],'b-')
        ax.plot([lag,lag],[fast-dfast,fast+dfast],'b-')
        ax.clabel(C,C.levels,inline=True,fmt ='%2.0f')
        # ax.set_title(r'Event {}. $\lambda$ 2 value = {}'.format(evt,lam2))




if __name__ == '__main__':
    print('This is discre.py')
    # p = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Accepted_SKS_SKKS_all.pairs'
    # path = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split'
    p = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Deng/Deng_events_SKS_SKKS.pairs'
    path = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Deng_events'
    t = Tester(p,path)
    lam2 =  t.pair_stack(p,path)
    p2 = write_lam2(p,lam2) # p2 contians lam2 values

    show_stacks(p2)
    #print(lam2)
    #plot_lam2(p.index.values,lam2)
