#! /usr/bin/env python
##############################
#   Program: synth.py
#
##############################
#   Author: J. Asplet
##############################
#   This program is designed to hold functions for reading SDBs containing synthetics
#   Make pairs file as required and to stack synthetics meausrements (i.e rather than
#   continue to hack away at sdb_analysis I wll treat synthetics as a special case And
#   writie them thier own special functions )
#
#
##############################
#   Import Statements
##############################
#   Standard Packages - all freely available
import os.path
import sys
import pandas as pd
import timeit
import shlex
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from stack import Stacker
###############################

class Synth:
    '''
    Class to hold the a Synthetics sdb. and then do some clever things with it ....
    '''
    def __init__(self,file,a=None):
        '''
        file - [str] The name of the synthetics SDB to read in
        '''
        date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        if file.split('.')[-1] == 'sdb':
            print('Sdb file read in')
            self.syn = pd.read_csv(file,delim_whitespace=True,converters=date_time_convert)
            self.nulls = self.syn[self.syn.Q <= -0.5] # df of clear nulls
            self.splits = self.syn[self.syn.Q >= 0.5] # df of clear splits
            self.unID = self.syn[(self.syn.Q > -0.5) & (self.syn.Q < 0.5)] # df of undeterminate events, in thoery there should be any ... ?
        elif file.split('.')[-1] == 'pairs':
            print('Pairs file read in. Can do plotting only. Also a is required')
            self.pairs =  pd.read_csv(file,delim_whitespace=True,converters=date_time_convert)
            self.a_ind = a
        self.spol= file.split('_')[0]  # [str] SP code for synthetics e.g. 'SP30'
        print(self.spol)
        self.F , self.T = np.meshgrid(np.arange(-90,95,5),np.arange(0,4.25,0.25))
        # Creates grids of fast and dt that we will need to plot the input synthetic splitting.
        self.lam2bar = [ ]
        self.lam2p1 = [ ]
        self.lam2p2 = [ ]
        nse = file.split('_')[2] # Pick out noise level from filename
        self.noise_lvl = ''.join(['Noise',nse])
        print(self.noise_lvl)

    def syn_in_v_out(self,save=False):
        '''Plot 2 subplot figure showing input synthetics (top) and the measurements made by sheba (bottom) '''
        fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True,figsize= (7,14))

        ax1.scatter(self.T,self.F,c='black',marker='.')
        ax1.set_xlim([0,4.0])
        ax1.set_ylim([-90,90])
        # ax1.set_xlabel(r'$\delta t$ (s)',fontsize=14)
        ax1.set_ylabel(r'$\Phi$ ( $\degree$)',fontsize=14 )
        ax1.set_title(r'$ \Phi$ v $\delta t$ space sampled by synthetics',fontsize=15)
        #Bottom Axis
        ax2.scatter(self.syn.TLAG,self.syn.FAST,marker='.',label='-0.7 < Q < 0.7')
        ax2.scatter(self.nulls.TLAG,self.nulls.FAST,marker='.',c='darkorange',label='Q < -0.7')
        ax2.set_xlim([0,4.0])
        ax2.set_ylim([-90,90])
        ax2.set_xlabel(r'$\delta t$ (s)',fontsize=16)
        ax2.set_ylabel(r'$\Phi$ ( $\degree$)' ,fontsize=16)
        ax2.set_title(r'Recovered $\Phi, \delta t$',fontsize=16)
        if save is True:
            plt.savefig('/Users/ja17375/Presentations/SYNTH_in_v_out.eps',format='eps',dpi=400)
        plt.tick_params(labelsize=14)
        plt.show()

    def grid_dSI(self,ax):
        ''' Plots dSI values (colurised) over the grid of fast, dt'''
        # fig = plt.figure() #figsize=(10,10))
        # ax = fig.add_subplot(111)
        f = self.F.ravel()
        l = self.T.ravel()
        dsi = self.pairs.D_SI.values.reshape(17,37)
        # C = ax.scatter(l,f,c=self.pairs.D_SI.values,marker='.',label='d_SI grid')
        # print(dsi.max())
        C = ax.contourf(self.T,self.F,dsi,18,levels=np.arange(0,dsi.max(),0.2),vmin=0.4,extend='max',cmap='viridis_r')
        # ax.contour(self.T,self.F,dsi,levels=[0.2],colors=['black'],linestyles='solid')
        ax.contour(self.T,self.F,dsi,levels=[0.4],colors=['black'],linestyles='solid')
        C.cmap.set_under('white')
        # Plot the singular A
        # ax.plot(self.a_lag,self.a_fast,'rx')
        ax.plot(l[self.a_ind],f[self.a_ind],'rx')
        # cbar1 = ax.colorbar(C,use_gridspec=False)
        # cbar1.set_label(r'$\Delta SI $',rotation=0)
        ax.set_ylabel(r'$\phi$ ($\degree$)',fontsize=14)
        ax.set_xlabel(r' $\delta t$ (s)',fontsize=14)
        ax.set_title(r'${}: \Delta SI$ for  $\delta t = {}, \phi = {}$ '.format(self.spol,l[self.a_ind],f[self.a_ind]))
        plt.tick_params(labelsize=14)
        return C
        # if save == True:
        #     print('dSI',self.spol,l[self.a_ind],f[self.a_ind])
        #     if f[self.a_ind] < 0:
        #         # print(abs(f[self.a_ind]))
        #         plt.savefig('/Users/ja17375/Presentations/{}_A_{:2.2f}_N{:03.0f}_dSI_grid.png'.format(self.spol,l[self.a_ind],abs(f[self.a_ind])),format='png',transparent=True,dpi=400)
        #         plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/Noise025/{}/{}_A_{:2.2f}_N{:03.0f}_dSI_grid.eps'.format(self.spol,self.spol,l[self.a_ind],abs(f[self.a_ind])),format='eps',transparent=True,dpi=400)
        #     else:
        #         plt.savefig('/Users/ja17375/Presentations/{}_A_{:2.2f}_{:03.0f}_dSI_grid.png'.format(self.spol,l[self.a_ind],f[self.a_ind]),format='png',transparent=True,dpi=400)
        #         plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/Noise025/{}/{}_A_{:2.2f}_{:03.0f}_dSI_grid.eps'.format(self.spol,self.spol,l[self.a_ind],f[self.a_ind]),format='eps',transparent=True,dpi=400)
        #     # plt.close('a')
        # plt.show()

    def grid_sigma2(self,ax):
        ''' Plots a grid of synthetic pairs coloured by whether they agree or disagree within 2-sigma'''
        fast = self.F.ravel() # Array of fast directions
        lag = self.T.ravel() # Array of lag times
        sig = np.zeros([17,37])

        for l in range(0,17): # We know there are 17 different lag times in our grid (0 - 4.0 at 0.25 spacing)
            #Define 2 sigma bound in fast direction for both phases
            lbt_P1 = self.pairs.TLAG_SKS[l] - 2*self.pairs.DTLAG_SKS[l]
            ubt_P1 = self.pairs.TLAG_SKS[l] + 2*self.pairs.DTLAG_SKS[l]
            lbt_P2 = self.pairs.TLAG_SKKS[l] - 2*self.pairs.DTLAG_SKKS[l]
            ubt_P2 = self.pairs.TLAG_SKKS[l] + 2*self.pairs.DTLAG_SKKS[l]

            for f in range(0,37): # We know there are 37 fast directions (-90 - 90 at 10 degree spacing)
                #Define 2 sigma bound in fast direction for both phases
                lbf_P1 = self.pairs.FAST_SKS[f] - 2*self.pairs.DFAST_SKS[f]
                ubf_P1 = self.pairs.FAST_SKS[f] + 2*self.pairs.DFAST_SKS[f]
                lbf_P2 = self.pairs.FAST_SKKS[f] - 2*self.pairs.DFAST_SKKS[f]
                ubf_P2 = self.pairs.FAST_SKKS[f] + 2*self.pairs.DFAST_SKKS[f]
                # Now test to see if the phases match or not
                if ((lbf_P2<= ubf_P1) & (lbf_P1 <= ubf_P2)) & ((lbt_P2 <= ubt_P1) & (lbt_P1 <= ubt_P2)):
                    # If test is passed that they match and the grid position is assigned a 1
                    sig[l,f] = 1
        # Plot grid of 1s and 0s. Potentially ned to find something better than contourf as it interpolats (?) slightly
        C = ax.contourf(self.T,self.F,sig)
        # C = ax.imshow(np.transpose(sig))
        ax.set_ylabel(r'$\phi$ ($\degree$)',fontsize=14)
        ax.set_xlabel(r' $\delta t$ (s)',fontsize=14)
        ax.plot(lag[self.a_ind],fast[self.a_ind],'rx')
        ax.set_title(r'{}: $2 \sigma$ for $\delta t = {}, \phi = {}$ '.format(self.spol,lag[self.a_ind],fast[self.a_ind]))
        # cbar1 = plt.colorbar(C,cax=ax,use_gridspec=True)
        return C

    def grid_lam2(self,ax):
        ''' Plots lam2 values (colurised) over the grid of fast, dt'''
        # fig = plt.figure()#figsize=(10,10))
        # ax = fig.add_subplot(111)
        f = self.F.ravel()
        l = self.T.ravel()
        lam2_bar = self.pairs.LAM2_BAR.values.reshape(17,37)
        lam2_p1 = self.pairs.LAM2_P1.values.min()
        lam2_p2 = self.pairs.LAM2_P2.values.min()
        l2_sum =  0.01 + (lam2_p1+lam2_p2)
        C = ax.contourf(self.T,self.F,lam2_bar,18,cmap='magma_r',extend='max')
        ax.contour(self.T,self.F,lam2_bar,levels=[l2_sum])
        # C.cmap.set_under('white')
        # Plot the singular A
        ax.plot(l[self.a_ind],f[self.a_ind],'rx')
        # cbar1 = ax.colorbar(C,use_gridspec=True)
        # cbar1.set_label(r'$\bar{\lambda_2} $',rotation=0)
        ax.set_ylabel(r'$\phi$ ($\degree$)',fontsize=14)
        ax.set_xlabel(r' $\delta t$ (s)',fontsize=14)
        ax.set_title(r'{}: $\lambda_2$ for $\delta t = {}, \phi = {}$ '.format(self.spol,l[self.a_ind],f[self.a_ind]))
        plt.tick_params(labelsize=14)
        return C
        # if save == True:
        #     print('Lam2', self.spol,l[self.a_ind],f[self.a_ind])
        #     if f[self.a_ind] < 0:
        #         # print(abs(f[self.a_ind]))
        #         plt.savefig('/Users/ja17375/Presentations/{}_A_{:2.2f}_N{:03.0f}_L2_grid.png'.format(self.spol,l[self.a_ind],abs(f[self.a_ind])),format='png',transparent=True,dpi=400)
        #         plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/Noise025/{}/{}_A_{:2.2f}_N{:03.0f}_L2_grid.eps'.format(self.spol,self.spol,l[self.a_ind],abs(f[self.a_ind])),format='eps',transparent=True,dpi=400)
        #     else:
        #         plt.savefig('/Users/ja17375/Presentations/{}_A_{:2.2f}_{:03.0f}_L2_grid.eps'.format(self.spol,l[self.a_ind],f[self.a_ind]),format='png',transparent=True,dpi=400)
        #         plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/Noise025/{}/{}_A_{:2.2f}_{:03.0f}_L2_grid.eps'.format(self.spol,self.spol,l[self.a_ind],f[self.a_ind]),format='eps',transparent=True,dpi=400)
        plt.close('all')
        plt.show()

    def plot_grids(self,save=False):
        '''Make a 4 panel plot of sigma2, deltaSI, lambda2bar grids along with a plot on the synthetics results'''
        # First, initialise fiugre instance and axes
        fig = plt.figure(figsize=(12,12))
        gs = gridspec.GridSpec(2,2) # Use gridspec to control axes allocation
        ax1 = plt.subplot(gs[0,0]) # Top left, axes to plot synthetics results onto
        ax2 = plt.subplot(gs[0,1]) # Top right, axes to plot 2sigma grid onto
        ax3 = plt.subplot(gs[1,0]) # Bottom left, axes to plot Delta SI grid onto
        ax4 = plt.subplot(gs[1,1]) # Bottom right, axes to plot lambda2 grid onto
        # Plot Synthetics as measured by sheba
        ax1.scatter(self.syn.TLAG,self.syn.FAST,marker='.',label='-0.7 < Q < 0.7')
        ax1.scatter(self.nulls.TLAG,self.nulls.FAST,marker='.',c='darkorange',label='Q < -0.7')
        ax1.set_xlim([0,4.0])
        ax1.set_ylim([-90,90])
        ax1.set_xlabel(r'$\delta t$ (s)',fontsize=16)
        ax1.set_ylabel(r'$\Phi$ ( $\degree$)' ,fontsize=16)
        ax1.set_title(r'Recovered $\Phi, \delta t$',fontsize=16)
        f = self.F.ravel()
        l = self.T.ravel()
        ax1.plot(l[self.a_ind],f[self.a_ind],'rx')
        # Plot 2 sigma grid
        C_2sigma = self.grid_sigma2(ax2)
        # fig.colorbar(C_2sigma,cax=ax2)
        # Plot delta SI grid
        C_dSI = self.grid_dSI(ax3)
        #Plot lamdba 2
        C_lam2 = self.grid_lam2(ax4)
        ########
        # Show plot
        if save == True:
            if f[self.a_ind] < 0:
                # print(abs(f[self.a_ind]))
                # plt.savefig('/Users/ja17375/Presentations/{}_A_{:2.2f}_N{:03.0f}_L2_grid.png'.format(self.spol,l[self.a_ind],abs(f[self.a_ind])),format='png',transparent=True,dpi=400)
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}/{}/{}_A_{:2.2f}_N{:03.0f}_4panel.png'.format(self.noise_lvl,self.spol,self.spol,l[self.a_ind],abs(f[self.a_ind])),format='eps',transparent=True,dpi=400)
            else:
                # plt.savefig('/Users/ja17375/Presentations/{}_A_{:2.2f}_{:03.0f}_L2_grid.eps'.format(self.spol,l[self.a_ind],f[self.a_ind]),format='png',transparent=True,dpi=400)
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}/{}/{}_A_{:2.2f}_{:03.0f}_4panel.png'.format(self.noise_lvl,self.spol,self.spol,l[self.a_ind],f[self.a_ind]),format='eps',transparent=True,dpi=400)

        plt.close('all')
        # plt.show()

    def plot_lamR(self,file,save=False):
        '''
        Plots a given raw lamR surface, to allow for comparison to the stacked synthetics grids (so we can see if the shape of the gridded lambda2 bar surfaces match
        the grid search output from SHEBA).
        '''
        lamR = np.loadtxt('/Users/Shear_Wave_Splitting/Sheba/Results/SYNTH/Stacks/{}'.format(file))
        L,F = np.meshgrid(np.linspace(0,4,num=lamR.shape[0]),np.linspace(-90,90,num=lamR.shape[1]))
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        C = ax.contourf(L,F,s30,20)
        ax.plot(2,-15,'rx')
        cbar = fig.colorbar(C)
        cbar.set_label(r'$\lambda _2$',rotation=0)
        ax.set_title(r'lamR surface for $\delta t = 2.0, \phi = -15$')
        ax.set_ylabel(r'$\phi (\degree) $ ')
        ax.set_xlabel(r'$\delta t$ (s)')
        if save == True:
            out = input('Input output file name (extension will be added) >')
            plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}.eps'.format(out),format='eps',transparent=True,dpi=400)
            print('Saved to /Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}.eps'.format(out))
        plt.show()

    def plot_in_v_out_joined(self,save=False):
        '''
        Plots the input grid of Fast and TLAG values along with the results from sheba. Pairs of points are joined by a line
        syn - df of sheba results
        a - start of range you want to plot (n.b. 0-18 give first column of dt )
        b - end of range. Defaults to all (i.e we set b to be the end of the synthetics file)
        '''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        f = self.F.ravel()
        l = self.T.ravel()
        # Now Plot pairs of events
        for i in range(0,len(self.syn)):
         x = [ l[i],self.syn.TLAG.values[i] ]
         y = [ f[i],self.syn.FAST.values[i]  ]
         if self.syn.Q.values[i] < -0.7:
             #pass
             ax.plot(x,y,'b-',color='darkorange',linewidth=0.5)
             ax.plot(l[i],f[i],'b.',color='darkorange',markersize=5)
             ax.plot(self.syn.TLAG.values[i],self.syn.FAST.values[i],'*',color='darkorange',markersize=5)
         else:
             ax.plot(x,y,'-',color='deepskyblue',linewidth=0.5)
             ax.plot(l[i],f[i],'r.',color='deepskyblue',markersize=5)
             ax.plot(self.syn.TLAG.values[i],self.syn.FAST.values[i],'*',color='deepskyblue',markersize=5)
             #pass
        ax.plot([-1,5],[30,30],'k--',linewidth=0.7)
        ax.plot([-1,5],[-60,-60],'k--',linewidth=0.7)
        ax.set_xlim([-0.5,4.5])
        ax.set_ylim([-95,95])
        ax.set_xlabel(r'$\delta t$ (s)')
        ax.set_ylabel(r'$\Phi $ ( $\degree$)')
        if save is True:
            fname = input('Enter output filename (eps extenstion will be added) > ')
            plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/{}.eps'.format(fname),format='eps',dpi=600)

        plt.show()


    def add_DSI(self):
        '''Calculate the difference in Splitting Intensity for each pair and add it to dataframe'''
        si_sks = self.pairs['SI(Pr)_x']
        si_skks = self.pairs['SI(Pr)_y']
        d_si = np.abs(si_sks-si_skks)
        self.pairs['D_SI'] = d_si
        #Delete SI cols as we dont need them any more ?
        # del self.pairs['INTENS_x']
        # del self.pairs['INTENS_y']


    def synth_stack(self,a,b):
        ''' Stack synthetics lamR surfaces for given a,b.
        a - [int] indicies of events you want to be the "SKS" (called SKS for sake of comparison to real data)
        b - [int] indicies of events you want to be the "SKKS" (called SKKS for sake of comparison to real data)
        noise_lvl [str] the Noise Level directory that the syntetics are contained it. (This is taken from self.noise_lvl)
        '''
    #Sanity Check
        if len(a) != len(b):
            Exception()
            #Throw excpetion if length of a,b dont match
        self.lam2bar = [ ]
        self.lam2p1 = [ ]
        self.lam2p2 = [ ]
        self.lam2sum = [ ]
        for i,k in enumerate(b):

            f1 = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/SYNTH/{}/{}/{}_3{:03d}001_120000_SYNTH.lamR'.format(self.noise_lvl,self.spol,self.spol,(a[i]+1)) # Add 1 to indecides becuase python goes from 0 628
            f2 = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/SYNTH/{}/{}/{}_3{:03d}001_120000_SYNTH.lamR'.format(self.noise_lvl,self.spol,self.spol,(b[i]+1)) # Whilst in the naing from BASH it goes from 1 to 629
            print(f1)
            print(f2)
            out = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/Stacks'.format(self.spol)
            outfile = '{}_3{:03d}{:03d}'.format(self.spol,a[i],b[i])
            Stk = Stacker(f1,f2,out,outfile)
            self.lam2bar.append(Stk.lam2_bar)
            self.lam2p1.append(Stk.lam2_sks)
            self.lam2p2.append(Stk.lam2_skks)
            # self.lam2sum = self.lam2p1 + self.lam2p2
        # Now we've done the stacks, add Lam2 to Synth pair_stack
        print('{} {} {}'.format(len(self.lam2bar),len(self.lam2p1),len(self.lam2p2))) #,len(self.lam2sum)))

        l2df = {'LAM2_BAR' : self.lam2bar, 'LAM2_P1' : self.lam2p1, 'LAM2_P2' : self.lam2p2} #,'LAM2_SUM' : self.lam2sum}
        ldf = pd.DataFrame(l2df)
        self.pairs[['LAM2_BAR','LAM2_P1','LAM2_P2',]] = ldf # 'LAM2_SUM' has been removed for now

    def synth_pairs(self,a,b,one_a=True,save=False,):
        '''makes synthetics pairs file sbased on user input (i,e you need to say which synthetics need pairing)
        a - [int] indicies of events you want to be the "SKS" (called SKS for sake of comparison to real data)
        b - [int] indicies of events you want to be the "SKKS" (called SKKS for sake of comparison to real data)
        noise_lvl [str] the Noise Level directory that the syntetics are contained it. This is the easiest way to get this information to make the pairs
        '''
        self.pairs = []
        # First lets select the relevent rows
        key = np.arange(0,len(a)) # Create a key so that we can merge rows without getting duplications
        # Slice out the rows we need to make the pairs
        A = self.syn.iloc[a].copy()
        B = self.syn.iloc[b].copy()
        # Resent the indexes so that they dont cuase problems in the merge
        A=A.reset_index(drop=True)
        B=B.reset_index(drop=True)
        #
        if one_a == True: # i.e we are stacking one point with the rest of the grid
            self.a_fast = A.FAST.values[0]
            self.a_lag = A.TLAG.values[0]
            self.a_ind = a[0]
            print(self.a_fast, self.a_lag)
        # Add the key
        A['key'] = pd.Series(key)
        B['key'] = pd.Series(key)
        # Perform the merge
        self.pairs = pd.merge(A,B,on=['key','TIME','STAT','STLA','STLO','EVLA','EVLO','EVDP','DIST','AZI','BAZ'],how='inner')
        #Relabel dulplicate columns to something meaningful
        relabel = {'FAST_x':'FAST_SKS', 'DFAST_x': 'DFAST_SKS','TLAG_x':'TLAG_SKS','DTLAG_x':'DTLAG_SKS','SPOL_x':'SPOL_SKS','DSPOL_x':'DSPOL_SKS',
              'WBEG_x':'WBEG_SKS','WEND_x':'WEND_SKS','EIGORIG_x':'EIGORIG_SKS','EIGCORR_x':'EIGCORR_SKS','Q_x':'Q_SKS','SNR_x':'SNR_SKS','NDF_x':'NDF_SKS',
              'FAST_y':'FAST_SKKS', 'DFAST_y': 'DFAST_SKKS','TLAG_y':'TLAG_SKKS','DTLAG_y':'DTLAG_SKKS','SPOL_y':'SPOL_SKKS','DSPOL_y':'DSPOL_SKKS',
              'WBEG_y':'WBEG_SKKS','WEND_y':'WEND_SKKS','EIGORIG_y':'EIGORIG_SKKS','EIGCORR_y':'EIGCORR_SKKS','Q_y':'Q_SKKS','SNR_y':'SNR_SKKS','NDF_y':'NDF_SKKS'}
        self.pairs.rename(relabel,axis='columns',inplace=True)
        # Create new DATE array. DATE for SYnthetics is of the form 3xxxyyy (mod of the yyyyjjj format where the leading 3 makes it obvious these arent real years)
        # xxx - index of the SKS synthetic
        # yyy - index of the SKKS synthetic
        D = ['3{:03d}{:03d}'.format(a[i],b[i]) for i in range(0,len(A))]
        # Add new DATEs
        self.pairs['DATE'] = pd.Series(D)
        self.pairs.drop(columns=['DATE_x','DATE_y','key'],inplace=True)
        # Get the DATE columns to the front of the array (where we would expect it normally)
        cols = self.pairs.columns.tolist()
        cols = cols[-1:] + cols[:-1] # Re-oreder columns so the new date is at the front
        self.pairs = self.pairs[cols]
        #Add d_SI values
        self.add_DSI()

        self.synth_stack(a,b)

        if save is True:
            f = self.F.ravel()
            l = self.T.ravel()
            # self.grid_dSI(save=True)
            # self.grid_lam2(save=True)
            plt.close('all')
            if f[self.a_ind] < 0:
                self.pairs.to_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/{}/{}/{}_A_{:2.2f}_N{:03.0f}.pairs'.format(self.noise_lvl,self.spol,self.spol,self.T.ravel()[self.a_ind],abs(self.F.ravel()[self.a_ind])),sep=' ')
            else:
                self.pairs.to_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/{}/{}/{}_A_{:2.2f}_{:03.0f}.pairs'.format(self.noise_lvl,self.spol,self.spol,self.T.ravel()[self.a_ind],self.F.ravel()[self.a_ind]),sep=' ')
        else:
            pass
            # self.grid_dSI()
            # self.grid_lam2()


    def l2_v_dSI(self):
        '''Plot lambda2 (bar) agaisnt splitting intensity'''
        fig,ax = plt.subplots(1,1,figsize=(6,6))
        ax.plot(self.pairs.LAM2_BAR,self.pairs.D_SI,'.')
        ax.set_xlabel(r'$\bar{\lambda_2} $ value', fontsize=14)
        ax.set_ylabel(r'$\Delta SI$ value',fontsize=14)
        plt.show()

    def l2bar_v_l2sum(self):
        ''' Plots lambda2 bar against lambda2 sum (ideally for a pair grid)'''
        fig,ax = plt.subplots(1,1,figsize=(7,7))
        mod = np.linspace(0,0.1,10)
        ax.plot(self.pairs.LAM2_SUM,self.pairs.LAM2_BAR,'k.')
        ax.plot(mod,mod*1.15)
        ax.set_xlabel(r'$\lambda_2^{P1} + \lambda_2^{P2}$')
        ax.set_ylabel(r'$\bar{\lambda_2}$')
        ax.set_xlim([0,0.005])
        ax.set_ylim([0,0.05])
        plt.show()

    def plot_snr(self):
        '''Plot Histograms of SNR for synthetics'''
        fig,ax = plt.subplots(1,1,figsize=(6,6))
        h = ax.hist(self.syn.SNR,20)
        ax.set_xlabel('Singal-to-Noise Ratio')
        ax.set_ylabel('Count')
        plt.show()


############ SYNTH CLASS ENDS #####################

def synth_l2_v_snr():

    pairs_p1_p1 = Synth('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/SP30_SNR_range_test_P1_P1.pairs').pairs
    pairs_p1_p2 = Synth('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/SP30_SNR_range_test_P1_P2.pairs').pairs
    pairs_p1_p3 = Synth('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/SP30_SNR_range_test_P1_P3.pairs').pairs
    pairs_p2_p2 = Synth('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/SP30_SNR_range_test_P2_P2.pairs').pairs
    pairs_p2_p3 = Synth('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/SP30_SNR_range_test_P2_P3.pairs').pairs
    pairs_p2_match = Synth('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/SP30_SNR_range_test_P2_match.pairs').pairs
    fig,(ax1,ax3) = plt.subplots(1,2,figsize= (12,6))
    #for pairs in pairs_list:

    snr_mod = np.arange(0.5,200,0.5)

    ax1.scatter(pairs_p1_p1.SNR,pairs_p1_p1.LAM2_BAR,marker='.',c='xkcd:blue',label=r'P1 - P1 pair "perfect" match')
    ax1.scatter(pairs_p2_match.SNR,pairs_p2_match.LAM2_BAR,marker='.',c='black',label=r'P2 - P2 pair "realistic" (0.2 diff in $\delta t$) match')
    ax1.scatter(pairs_p2_p2.SNR,pairs_p2_p2.LAM2_BAR,marker='.',c='xkcd:green',label=r'P2 - P2 pair "perfect" match')
    ax1.scatter(pairs_p1_p2.SNR,pairs_p1_p2.LAM2_BAR,marker='.',c='xkcd:goldenrod',label=r'P1 - P2 pair (discrepant)')
    #ax1.scatter(pairs_p1_p3.SNR,pairs_p1_p3.LAM2_BAR,marker='.',c='xkcd:pink',label=r'P1 - P3 pair (discrepant)')
    ax1.scatter(pairs_p2_p3.SNR,pairs_p2_p3.LAM2_BAR,marker='.',c='xkcd:red',label='P2 - P3 pair')
    ax1.set_title(r'$\bar{\lambda_2}$ for split pairs against SNR')
    ax1.set_xlabel('SNR of SYN1')
    ax1.legend()
    # Make a log-log plot for lam2 v SNR
    ax3.loglog(pairs_p1_p1.SNR,pairs_p1_p1.LAM2_BAR,'.',color='xkcd:blue',label=None)
    ax3.loglog(pairs_p2_p2.SNR,pairs_p2_p2.LAM2_BAR,'.',color='xkcd:green',label=None)
    ax3.loglog(pairs_p2_match.SNR,pairs_p2_match.LAM2_BAR,'.',color='black',label=None)
    ax3.loglog(pairs_p1_p2.SNR,pairs_p1_p2.LAM2_BAR,'.',color='xkcd:goldenrod',label=None)
    #ax3.loglog(pairs_p1_p3.SNR,pairs_p1_p3.LAM2_BAR,'.',color='xkcd:pink',label=None)
    ax3.loglog(pairs_p2_p3.SNR,pairs_p2_p3.LAM2_BAR,'.',color='xkcd:red',label=None)
    # Fit a straight line to THIS (log-log) data
    m2, c2 = np.polyfit(np.log(pairs_p1_p1.SNR),np.log(pairs_p1_p1.LAM2_BAR),1)
    # Calculate line
    y_fit = m2*np.log(snr_mod) + c2
    #print(max(SNR),min(SNR))
    # ax3.plot(snr_mod,10**y_fit,'k--')
    ax3.set_title(r'loglog plot of $\bar{\lambda_2}$ for split pairs against SNR')
    # Plot models on ax1
    print(m2,c2)
    ax3.text(5,0.1,r'$log(\lambda _2$) = {:4.2f} * log(SNR) + {:4.2f}'.format(m2,c2),transform=ax3.transAxes)
    ax3.set_ylabel(r'log($\lambda _2$)')
    ax3.set_xlabel('SNR')

    y_mod = (m2 * np.log(snr_mod)) + c2

    # Set some limits

    #f1, = ax1.plot(snr_mod,y_fit1,'--',label='Semilog fit')
    ax1.plot(snr_mod,np.exp(y_mod),'k--',label='Model')
    ax1.plot(snr_mod,np.exp(y_mod)+0.003, 'k-.')
    f1, = ax3.loglog((snr_mod),np.exp(y_mod), 'k--',label='Model for "Perfect" matches')
    f2, = ax3.loglog(snr_mod,np.exp(y_mod)+0.003 ,'k-.',label=r'$\bar{\lambda _2} = A. SNR^{-b} + \kappa_{\bar{\lambda_2}}$ Model')
    #ax1.plot([10,10],[0,L2.max()],'-k')
    ax1.set_ylim([-0.01,0.3])
    ax1.set_xlim([1,55])
    ax3.set_xlim([1,55])
    ax3.legend([f1,f2],['Model for Perfect Matches',r'$\bar{\lambda _2} = A. SNR^{-b} + \kappa_{\bar{\lambda_2}}$ Model'])
    #f4 = ax3.plot(snr_mod,y_mod,'k-')
    #plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/LAM2_v_SNR_w_models.eps',format='eps',dpi=400)
    #plt.savefig('LAM2_v_SNR_matches.eps',format='eps',dpi=400)
    plt.show()






def synth_l2bar_v_l2(self):
    '''
    Form Ipython, needs recaasting into Class stcutre
    '''
    fig,ax = plt.subplots(1,1,figsize=(7,7))
    C1 = ax.scatter((pairs_p1_p1.LAM2_P1+pairs_p1_p1.LAM2_P2),pairs_p1_p1.LAM2_BAR,marker='.',c='black',label=r'S1 - S1 pair "perfect" match',vmax=20.0)
    ax.scatter((pairs_p2_match.LAM2_P1+pairs_p2_match.LAM2_P2),pairs_p2_match.LAM2_BAR,marker='.',c='green',label=r'S2 - S2 pair "realistic" (0.2 diff in $\delta t$) match')
    #ax.scatter((pairs_p2_p2.LAM2_P1+pairs_p2_p2.LAM2_P2),pairs_p2_p2.LAM2_BAR,marker='.',c='xkcd:green',label=r'P2 - P2 pair "perfect" match')
    C2 = ax.scatter((pairs_p1_p2.LAM2_P1+pairs_p1_p2.LAM2_P2),pairs_p1_p2.LAM2_BAR,marker='.',c='red',label=r'S1 - S2 pair (discrepant)',vmax=20.0)
    #ax1.scatter(pairs_p1_p3.SNR,pairs_p1_p3.LAM2,marker='.',c='xkcd:pink',label=r'P1 - P3 pair (discrepant)')
    C3 = ax.scatter((pairs_p2_p3.LAM2_P1+pairs_p2_p3.LAM2_P2),pairs_p2_p3.LAM2_BAR,marker='.',c='blue',label='S2 - S3 pair',vmax=20.0)
    ax.set_xlabel(r'$\lambda_2^{P1} + \lambda_2^{P2}$')
    ax.set_ylabel(r'$\bar{\lambda_2}$')
    s1 = r'S1: $\phi = -10, \delta t = 2.0$'
    s2 =r'S2: $\phi = 70, \delta t = 1.5 $'
    s3 = r'S3: $\phi = 50, \delta t = 2.5 $'
    s2m = r'S2_match: $ \phi = 70, \delta t = 1.3$.'
    mod = np.linspace(0,0.8,10)
    ax.plot(mod,mod,'k--',label=r'$\bar{\lambda_2} = \lambda_2^{P1} + \lambda_2^{P2}$')
    ax.plot(mod,mod+0.1,'k-.',label=r'$\bar{\lambda_2} = (\lambda_2^{P1} + \lambda_2^{P2}) + 0.1$')
    ax.legend()
    ax.text(0.125,0.1,s1,fontsize=10,wrap=True)
    ax.text(0.125,0.092,s2,fontsize=10,wrap=True)
    ax.text(0.125,0.084,s2m,fontsize=10,wrap=True)
    ax.text(0.125,0.076,s3,fontsize=10,wrap=True)
    #ax.text(0.2,0.08,s_cont,fontsize=10)
    ax.set_ylim([0,0.2])
    ax.set_xlim([0,0.2])
    # Add label detailing the splitting parmeters
    plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/Lam2bar_v_lam2sum_synthetics.eps',format='eps',dpi=400)
    plt.show()

def mk_syn_pairs(a,b):
     '''
     a = sdb for the first phase to be stacked
     b = second sdb
     '''
     pairs = pd.merge(a,b,on=['DATE','TIME','STAT','STLA','STLO','EVLA','EVLO','EVDP','DIST','AZI','BAZ'],how='inner')
     relabel = {'FAST_x':'FAST_P1', 'DFAST_x': 'DFAST_P1','TLAG_x':'TLAG_P1','DTLAG_x':'DTLAG_P1','SPOL_x':'SPOL_P1','DSPOL_x':'DSPOL_P1',
               'WBEG_x':'WBEG_P1','WEND_x':'WEND_P1','EIGORIG_x':'EIGORIG_P1','EIGCORR_x':'EIGCORR_P1','Q_x':'Q_P1','SNR_x':'SNR_P1','NDF_x':'NDF_P1',
               'FAST_y':'FAST_P2', 'DFAST_y': 'DFAST_P2','TLAG_y':'TLAG_P2','DTLAG_y':'DTLAG_P2','SPOL_y':'SPOL_P2','DSPOL_y':'DSPOL_P2',
               'WBEG_y':'WBEG_P2','WEND_y':'WEND_P2','EIGORIG_y':'EIGORIG_P2','EIGCORR_y':'EIGCORR_P2','Q_y':'Q_P2','SNR_y':'SNR_P2','NDF_y':'NDF_P2'}
     pairs.rename(relabel,axis='columns',inplace=True)
     pairs['SNR'] = pairs.SNR_P1 + pairs.SNR_P2
     pairs['D_SI'] = np.abs(pairs['SI(Pr)'] - pairs['SI(Pr)']) # Only use projection for now
     # del pairs['INTENS_x']
     # del pairs['INTENS_y']
     ldf = synth_stack(pairs.DATE.values)
     #print(ldf)
     pairs[['LAM2_BAR','LAM2_P1','LAM2_P2']] = ldf
     return pairs


def synth_stack(a):
     ''' Stack  lamR surfaces for a set of synthetics of equal length'''
     lam2_bar = [ ]
     lam2_P1 = [ ]
     lam2_P2 = [ ]
     for i,k in enumerate(a):
         f1 = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/SYNTH/SNR_test/P3/SP30_{}_120000_SYNTH.lamR'.format(a[i]) # Select synth with correct date from direcotry P1
         f2 = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/SYNTH/SNR_test/P3/SP30_{}_120000_SYNTH.lamR'.format(a[i])  # Select synth from P2
         print(f1)
         print(f2)

         out = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/SNR_test/Stacks'
         outfile = 'SP30_{:07d}'.format(a[i])
         Stk = Stacker(f1,f2,out,outfile)
         lam2_bar.append(Stk.lam2_bar)
         lam2_P1.append(Stk.lam2_sks)
         lam2_P2.append(Stk.lam2_skks)
     # Now we've done the stacks, add Lam2 to Synth pair_stack
     l2df = {'LAM2_BAR' : lam2_bar,'LAM2_P1' : lam2_P1,'LAM2_P2' : lam2_P2 }
     ldf = pd.DataFrame(l2df)
     return ldf
