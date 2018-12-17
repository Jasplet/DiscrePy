#! /anaconda3/envs/splitwavepy/bin/python
##########################################
# stack.py
##########################################
# Program that stacks lam2 surfaces for SKS SKKS pairs using
# the sheba_stack fortran program. This is used to look for
# discrpeent splitting
##########################################
import obspy as ob
import numpy as np
import pandas as pd
import subprocess as sub
import os
from os import path
import shutil
import matplotlib.pyplot as plt
##########################################

# Identify pairs that we want to stack
# or have them passed in as arguements?

# cd to sheba runs directory
# cp lam2 files to "right place"
# generate sheba_stack.in
# subprocess call sheba stack,

if __name__ == '__main__':
    '''
    Main level, expects 2 arguements for the lam2 surfaces in questions
    '''


class Stacker:
    '''
    Class for doing the stacking
    '''

    def __init__(self,lam2_sks,lam2_skks,outpath,outfile=None,stk_type='man'):#,syn='False',i=None):
        '''
        Initialises class, checks if lam2 surfaces exist
        lam2_sks - [str] full path to lam2 surface for SKS
        lam2_skk - [str] full path to lam2 surface for SKKS
        fstem - [str]
        These should be provided as FULL PATHS !
        '''
        # print('Starting Stacker')

        if os.path.isfile(lam2_sks) is False:
            raise NameError('Lambda 2 (for sks) provided does not exist')
        elif os.path.isfile(lam2_skks) is False:
            raise NameError('Lambda 2 (for skks) provided does not exist')
        # else:
            # print('Lambda 2 surfaces exist')
    #   make input paths attribute, we need these full paths for man stacking mode
        self.sks_path = lam2_sks
        self.skks_path = lam2_skks
        #   Isolate filestems of lambda 2 surfaces (for SHEBA mode)
        self.sks = lam2_sks.split('/')[-1]
        self.skks = lam2_skks.split('/')[-1]
    #   make output filestem (but cutting off the phase extension from self.sks)

        if outfile == None:
            self.out = '_'.join(self.sks.split('_')[:-1])
            # print('Outfile:', self.out)
        else:
            print(outfile)
            self.out = outfile

        #Copy .lamR files to our Stacks directory so we still have them for plotting if I decide to Purge the Runs directory
        # print('Outpath:',outpath)
        self.copy_files(lam2_sks,lam2_skks,outpath)

        self.outfile = '{}/{}'.format(outpath,self.out)
    #   make arrays of dt and fast that we can use to identify solution of the stack (for man mode)
        self.T = np.arange(0,4.025,0.025)
        self.F = np.arange(-90,91,1)

#       Perform stack
        if stk_type == 'sheba':

            path_stem = lam2_sks.split('/')[0:7] #+ [fstem]
            # print(path_stem)
            self.path = '/'.join(path_stem)
            if os.path.isdir(self.path) is False:
                os.mkdir(self.path)
            #Copy lam2 files to where we want to work on them
            self.copy_files(lam2_sks,lam2_skks,self.path)
    #       Make infile
            self.make_infile()
            self.stack_sheba()
            self.collect()
        elif stk_type == 'man':
            self.stack_manual()



    def stack_sheba(self):
        print('Stacking')
        p=sub.Popen(['sheba_stack'],stdout = sub.PIPE,
                                    stdin  = sub.PIPE,
                                    stderr = sub.STDOUT,
                                    cwd = self.path,
                                    encoding='utf8')

        p.communicate('-wgt one')

    def stack_manual(self):
        # print('Manual stacker, function under construction, this sis a placeholder command')

        self.sks_lamR = np.loadtxt(self.sks_path)
        self.skks_lamR = np.loadtxt(self.skks_path)

        #perform stack by adding surfaces together.
        # No weighting applied
        self.stk = (self.sks_lamR + self.skks_lamR)

        # find min lam2 value - for both phases and for the stacked surface
        self.lam2_sks = self.sks_lamR.min()
        self.lam2_skks = self.skks_lamR.min()
        self.lam2_bar = self.stk.min()
        # find its location

        jf,jt  = np.unravel_index(self.stk.argmin(),self.stk.shape)
        print('Min Lam2 of stack is {}, located at dt = {}  and phi = {}'.format(self.lam2_bar,self.T[jt],self.F[jf]))
        # print(self.outfile)
        # write out stack
        # print(self.stk.shape)

        np.savetxt('{}.lamSTK'.format(self.outfile),self.stk,fmt='%.5f')

        self.sol = self.lam2_bar


    def collect(self):
        ''' Collects lambda2 value and solution of the stacked surface '''
        with open('{}/sheba_stack.sol'.format(self.path)) as reader:
            header = reader.readline()
            sol = reader.readline()
            f = sol.strip('\n').split(' ')[0]
            df = sol.strip('\n').split(' ')[1]
            l = sol.strip('\n').split(' ')[2]
            dl = sol.strip('\n').split(' ')[3]
            lam2 = sol.strip('\n').split(' ')[-1]
            # print('lambda 2 value is {}'.format(lam2))
            self.sol = [f,df,l,dl,float(lam2)]

    def copy_files(self,sks,skks,path):
        ''' Copies lambda 2 files to corrrect place for stacking '''
        sub.call(['cp',sks,path])
        sub.call(['cp',skks,path])


    def make_infile(self):
        ''' Makes the sheba_stack.in file thats required '''
        with open('{}/sheba_stack.in'.format(self.path),'w') as writer:
            writer.write('{} \n'.format(self.sks))
            writer.write(self.skks)

        # print('sheba_stack.in written to {}'.format(self.path))


########

def plot_stack(paths):
    ''' Function to read a sheba stack .sol and. err file and plot the stacked SKS and SKKS surfaces -  Adapted from plot_sheba_stack.m by J Wookey'''

    for i,path in enumerate(paths):
        p = path.split('/')
        print(p)
    # Read solution
        with open('{}/sheba_stack.sol'.format(path),'r') as reader:
            head = reader.readline()  #Reads headers
            S = reader.readline().split() # Reads solution

            fast,dfast = float(S[0]), float(S[1])
            lag,dlag = float(S[2]),float(S[3])
            nsurf = float(S[4])
            lag_step = float(S[5])
            lam2 = S[6]
            print(lam2)
    # Read surface
    err = np.loadtxt('{}/sheba_stack.err'.format(path))


    nfast,nlag = err.shape ;

    lag_max = (nlag) * lag_step ;
    [T,F] = np.meshgrid(np.arange(0,lag_max,lag_step),np.arange(-90,91,1)) ;
    fig = plt.figure(i)
    C = plt.contour(T,F,err,[0.5,1,2,3,4,5,10,15,20,25,30,40,50],colors='k')
    plt.ylabel(r'Fast,$\phi$, (deg)')
    plt.xlabel(r'Lag ,$\delta$ t, (sec)')
    plt.plot([lag-dlag,lag+dlag],[fast,fast],'b-')
    plt.plot([lag,lag],[fast-dfast,fast+dfast],'b-')
    plt.clabel(C,C.levels,inline=True,fmt ='%2.0f')
    plt.title(r'Event {}. $\lambda$ 2 value = {}'.format(p[8],lam2))

    plt.show()
# %   figure
#    [C,h] = contour(TLAG,FAST,ERR,[1,1],'k-') ;
#    set(h,'LineWidth',2)
#    hold on
#    [C,h] = contour(TLAG,FAST,ERR,[1,2,3,4,5,10,15,20,50,100],'k-') ;
#    clabel(C,h) ;
#
#    plot([tlag-dtlag;tlag+dtlag],[fast;fast],'b-')
#    plot([tlag ;tlag ],[fast-dfast;fast+dfast],'b-')
#
#    axis([0 tlag_max -90 90]) ;
#    xlabel('TLAG (SEC)','FontSize',12,'FontWeight','bold') ;
#    ylabel('FAST (DEG)','FontSize',12,'FontWeight','bold') ;
#
#
#    figure
#
#    pcolor(TLAG,FAST,ERR)
#    colorbar
