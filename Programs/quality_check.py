#! /usr/bin/env python
######################
# quality_check.py
######################
# A program to loop through and inspect SWS results, plot the QA result from
# Splitwavepy and print the measured SKS and SKKS splitting from Sheba
# A measure of measure of data quality and discernable discrepancy between results is then assigned
# to each event. This will allow bad results to by removed
#####################################################################################################
# Some Imports
import pandas as pd
import splitwavepy as sw
import obspy
import matplotlib.pyplot as plt
from summary_plot import plotall
import os
import SKS_SKKS_qa


class Inspecter:

    def __init__(self, pairs,dir):
        self.pairs =  pd.read_csv('{}/{}'.format(dir,sys.argv[1]),delim_whitespace=True,converters={'TIME': lambda x: str(x)})
        self.result_path = dir
        # self.out = 


    def loop_thru_pairs(self):
        '''Loops through the pairs and checks if there is a set of QA results'''
        for i,row in self.pairs.iterrows():
            filestem = '{}_{}_{}'.format(self.pairs.STAT[i],self.pairs.DATE[i],self.pairs.TIME[i])



if __name__ == '__main__':
    print('Hello World, this is quality_check.py. You are running me from the command line')

    pair_file = sys.argv[1]
    results_dir = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split'

    if len(sys.argv) is 1:
        print('No inputs detected. Please input them now')


    Review = Inspecter(pair_file,results_dir)
