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
        self.lam2 = [ ]


    def syn_in_v_out(self):
        '''Plot 2 subplot figure showing input synthetics (top) and the measurements made by sheba (bottom) '''
        fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,sharex=True,figsize= (8,8))

        ax1.scatter(self.T,self.F,c='black',marker='.')
        ax1.set_xlim([0,4.0])
        ax1.set_ylim([-90,90])
        ax1.set_xlabel(r'$\delta t$ (s)')
        ax1.set_ylabel(r'$\Phi$ ( $\degree$)' )
        ax1.set_title(r'$ \Phi$ v $\delta t$ space sampled by synthetics')
        #Bottom Axis
        ax2.scatter(self.syn.TLAG,self.syn.FAST,marker='.',label='-0.7 < Q < 0.7')
        ax2.scatter(self.nulls.TLAG,self.nulls.FAST,marker='.',c='darkorange',label='Q < -0.7')
        ax2.set_xlim([0,4.0])
        ax2.set_ylim([-90,90])
        ax2.set_xlabel(r'$\delta t$ (s)')
        ax2.set_ylabel(r'$\Phi$ ( $\degree$)' )
        ax2.set_title('Shear-wave splitting measured from Synthetics by Sheba')
        plt.show()

    def grid_dSI(self,save=False):
        ''' Plots dSI values (colurised) over the grid of fast, dt'''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        f = self.F.ravel()
        l = self.T.ravel()
        dsi = self.pairs.D_SI.values.reshape(17,37)
        # C = ax.scatter(l,f,c=self.pairs.D_SI.values,marker='.',label='d_SI grid')
        C = ax.contourf(self.T,self.F,dsi,18,vmin=0,extend='max')
        C.cmap.set_under('white')
        # Plot the singular A
        # ax.plot(self.a_lag,self.a_fast,'rx')
        ax.plot(l[self.a_ind],f[self.a_ind],'rx')
        cbar1 = fig.colorbar(C,use_gridspec=False)
        cbar1.set_label(r'$\delta SI $',rotation=270)
        ax.set_ylabel(r'$\phi$ ($\degree$)')
        ax.set_xlabel(r' $\delta t$ (s)')
        ax.set_title(r'${}: \Delta SI$ for  $\delta t = {}, \phi = {}$ '.format(self.spol,l[self.a_ind],f[self.a_ind]))
        if save == True:
            print('dSI',self.spol,l[self.a_ind],f[self.a_ind])
            if f[self.a_ind] < 0:
                # print(abs(f[self.a_ind]))
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}_A_{:2.2f}_N{:03.0f}_dSI_grid.eps'.format(self.spol,l[self.a_ind],abs(f[self.a_ind])),format='eps',transparent=True,dpi=400)
            else:
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}_A_{:2.2f}_{:03.0f}_dSI_grid.eps'.format(self.spol,l[self.a_ind],f[self.a_ind]),format='eps',transparent=True,dpi=400)
            plt.close('a')
        else:
            plt.show()

    def grid_lam2(self,save=False):
        ''' Plots dSI values (colurised) over the grid of fast, dt'''
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        f = self.F.ravel()
        l = self.T.ravel()
        lam2 = self.pairs.LAM2.values.reshape(17,37)
        C = ax.contourf(self.T,self.F,lam2,18,cmap='magma_r',vmin=0,vmax=0.1,extend='both')
        C.cmap.set_over('black')

        #C = ax.scatter(l,f,self.LAM2,marker='.',label='d_SI grid')
        # Plot the singular A
        # ax.plot(self.a_lag,self.a_fast,'rx')
        ax.plot(l[self.a_ind],f[self.a_ind],'rx')
        cbar1 = fig.colorbar(C,use_gridspec=True)
        cbar1.set_label(r'$\bar{\lambda_2} $',rotation=0)
        ax.set_ylabel(r'$\phi$ ($\degree$)')
        ax.set_xlabel(r' $\delta t$ (s)')
        ax.set_title(r'{}: $\lambda _2$ for $\delta t = {}, \phi = {}$ '.format(self.spol,l[self.a_ind],f[self.a_ind]))
        if save == True:
            print('Lam2', self.spol,l[self.a_ind],f[self.a_ind])
            if f[self.a_ind] < 0:
                # print(abs(f[self.a_ind]))
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}_A_{:2.2f}_N{:03.0f}_L2_grid.eps'.format(self.spol,l[self.a_ind],abs(f[self.a_ind])),format='eps',transparent=True,dpi=400)
            else:
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/SynthStacks/{}_A_{:2.2f}_{:03.0f}_L2_grid.eps'.format(self.spol,l[self.a_ind],f[self.a_ind]),format='eps',transparent=True,dpi=400)
            plt.close('a')
        else:
            plt.show()

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
        si_sks = self.pairs.INTENS_x
        si_skks = self.pairs.INTENS_y
        d_si = np.abs(si_sks-si_skks)
        self.pairs['D_SI'] = d_si
        #Delete SI cols as we dont need them any more ?
        del self.pairs['INTENS_x']
        del self.pairs['INTENS_y']


    def synth_stack(self,a,b):
        ''' Stack synthetics lamR surfaces for given a,b. '''
    #Sanity Check
        if len(a) != len(b):
            Exception()
            #Throw excpetion if length of a,b dont match

        for i,k in enumerate(b):

            f1 = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/SYNTH/{}/{}_3{:03d}001_120000_SYNTH.lamR'.format(self.spol,self.spol,(a[i]+1)) # Add 1 to indecides becuase python goes from 0 628
            f2 = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/SYNTH/{}/{}_3{:03d}001_120000_SYNTH.lamR'.format(self.spol,self.spol,(b[i]+1)) # Whilst in the naing from BASH it goes from 1 to 629
            print(f1)
            print(f2)
            out = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/Stacks'.format(self.spol)
            outfile = '{}_3{:03d}{:03d}'.format(self.spol,a[i],b[i])
            Stk = Stacker(f1,f2,out,outfile)
            self.lam2.append(Stk.lam2)

        # Now we've done the stacks, add Lam2 to Synth pair_stack
        l2df = {'LAM2' : self.lam2}
        ldf = pd.DataFrame(l2df)
        self.pairs['LAM2'] = ldf

    def synth_pairs(self,a,b,one_a=True,save=False):
        '''makes synthetics pairs file sbased on user input (i,e you need to say which synthetics need pairing)
        a - [int] indicies of events you want to be the "SKS" (called SKS for sake of comparison to real data)
        b - [int] indicies of events you want to be the "SKKS" (called SKKS for sake of comparison to real data)
        '''
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
            if f[self.a_ind] < 0:
                self.pairs.to_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/{}_lowSNR_A_{:2.2f}_N{:03.0f}_B_grid.pairs'.format(self.spol,self.T.ravel()[self.a_ind],abs(self.F.ravel()[self.a_ind])),sep=' ')
            else:
                self.pairs.to_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/{}__lowSNR_A_{:2.2f}_{:03.0f}_B_grid.pairs'.format(self.spol,self.T.ravel()[self.a_ind],self.F.ravel()[self.a_ind]),sep=' ')
