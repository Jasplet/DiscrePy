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

    def plot_baz():
        ''' Make plot of Fast and Lag v BAZ for SKS and SKKS'''

    def plot_dSI():
        '''plot a histogram dSI for the bin'''

    def plot_lam2():
        '''plot a histogram of LAM2 values for the bin'''

    def avg_lam2():
        '''return average lambda 2 value assuming a guassian distribution'''

    def avg_dSI():
        '''return average lambda 2 value assuming a gaussian distribution in the bin'''

    def avg_splitting():
        '''return average FAST,LAG for the bin. Assuming a gaussian distribution'''



if __name__ == '__main__':
  print('Hello I am bin_analysis.py! I am not currently deisgn ed to work from the command line. Please trying importing me into your script or using IPYTHON  :-) ')
