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

    def __init__(self,lam2_sks,lam2_skks):
        '''
        Initialises class, checks if lam2 surfaces exist

        These should be provided as FULL PATHS !
        '''

        if os.path.isfile(lam2_sks) is False:
            raise NameError('Lambda 2 (for sks) provided does not exist')
        elif os.path.isfile(lam2_skks) is False:
            raise NameError('Lambda 2 (for sks) provided does not exist')
        else:
            print('Lambda 2 sufraces exist')
    #   Isolate filestems of lambda 2 surfaces

        self.sks = lam2_sks.split('/')[-1]
        self.skks = lam2_skks.split('/')[-1]

        path_stem = lam2_sks.split('/')[0:8]
        self.path = '/'.join(path_stem)
        #Copy lam2 files to where we want to work on them
        self.copy_files(lam2_sks,lam2_skks)
#       Make infile
        self.make_infile()
#       Perform stack
        self.stack()
        self.collect()

    def stack(self):
        print('Stacking')
        p=sub.Popen(['sheba_stack'],stdout = sub.PIPE,
                                    stdin  = sub.PIPE,
                                    stderr = sub.STDOUT,
                                    cwd = self.path,
                                    encoding='utf8')

        p.communicate('-wgt one')

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
            print('lambda 2 value is {}'.format(lam2))
            self.sol = [f,df,l,dl,float(lam2)]

    def copy_files(self,sks,skks):
        ''' Copies lambda 2 files to corrrect place for stacking '''
        sub.call(['cp',sks,self.path])
        sub.call(['cp',skks,self.path])


    def make_infile(self):
        ''' Makes the sheba_stack.in file thats required '''
        with open('{}/sheba_stack.in'.format(self.path),'w') as writer:
            writer.write('{} \n'.format(self.sks))
            writer.write(self.skks)

        print('sheba_stack.in written to {}'.format(self.path))

########

def plot_stack():
    ''' Function to read a sheba stack .sol and. err file and plot the stacked SKS and SKKS surfaces -  Adapted from plot_sheba_stack.m by J Wookey'''

    # Read solution
    with open('sheba_stack.sol','r') as reader:
        head = reader.readline()  #Reads headers
        S = reader.readline().split() # Reads solution

        fast,dfast = S[0], S[1]
        lag,dlag = S[2],S[3]
        nsurf = S[4]
        lag_step = S[5]
        lam2 = S[6]

    # Read surface
    err = np.loadtxt('sheba_stack.err')


    nfast,nlag = err.shape ;

    lag_max = (nlag-1) * lag_step ;
    [T,F] = np.meshgrid(np.arange(0,lag_max,lag_step),np.arange(-90,90,1)) ;
#
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