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
import os.path
import SKS_SKKS_qa
import sys

class Inspecter:

    def __init__(self, pairs,dir):
        self.pairs =  pd.read_csv('{}/{}'.format(dir,pairs),delim_whitespace=True,converters={'TIME': lambda x: str(x)})
        self.result_path = dir
        self.qa_dir = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/SplitWavePy'
        self.accepted_i = []
        self.fs = []

    def write_qa(self,filestem, qual, disc,mode='a+'):
        ''' Writes the event filestem and the quality/discrepancy rating to a textfile'''

        with open('{}/Inspection_Results.txt'.format(self.result_path),mode) as outfile:
            outfile.write('{} {} {}\n'.format(filestem,qual,disc) )



    def run_inspection(self):
        if os.path.isfile('{}/Inspection_Results.txt'.format(self.result_path)):
            chk = input('Inspection_Results.txt already exists, [o]verwrite or [c]ontinue? \n >')
            if chk is 'c':
                with open('{}/Inspection_Results.txt'.format(self.result_path)) as reader:
                    head= reader.readline()
                    for line in reader.readlines():
                        self.fs = line.split(' ')[0]
                        self.qual = line.split(' ')[1]
                        self.disc = ' '.join(line.split(' ')[2:])

                self.loop_thru_pairs()
        else:
            print('Initialise Inspection_Results.txt')
            self.write_qa('Filestem', 'Data Quality', 'Discrepancy',mode='w+')
            self.loop_thru_pairs()

    def loop_thru_pairs(self):
        '''Loops through the pairs and checks if there is a set of QA results'''
        # Write HEaderline to outfile
        for i,row in self.pairs.iterrows():
            filestem = '{}_{}_{}'.format(row.STAT,row.DATE,row.TIME)
            if filestem in self.fs:
                print("Event already inspected. Qual {} Disc {}".format(self.qual[i],self.disc[i]))
                qual = self.qual[i]
                disc = self
            print('Event {} Date: {} Time: {}, Stat: {}'.format(i,row.DATE,row.TIME,row.STAT))

            if row.SNR_SKS >= 10 or row.SNR_SKKS >= 10:
                #Test to see if Signal-to-Noise is too high
                print('SNR too high, auto-reject')
                qaul = 'p'
                disc = 'SNR too high'
            else:

                sr = (row.FAST_SKS,row.DFAST_SKS,row.TLAG_SKS,row.DTLAG_SKS,row.FAST_SKKS,row.DFAST_SKKS,row.TLAG_SKKS,row.DTLAG_SKKS)
                print('SHEBA RESULTS: SKS: phi = {} +/- {} dt = {} +/- {}. SKKS: phi = {} +/- {}, dt = {} +/- {}'.format(sr[0],sr[1],sr[2],sr[3],sr[4],sr[5],sr[6],sr[7]))
                #Test if QA results (from SWP) exist
                if (os.path.isfile('{}/{}_sks.eigm'.format(self.qa_dir,filestem)) and os.path.isfile('{}/{}_skks.xcrm'.format(self.qa_dir,filestem))) is False:
                    print('QA results do not exist, generating')
                    SKS_SKKS_qa.measure_sks_skks(filestem,self.qa_dir,[row.WBEG_SKS,row.WEND_SKS,row.WBEG_SKKS,row.WEND_SKKS])
                    # Now plot the TransM,Xcross and Eigm results
                    # fig = plt.figure(figsize=(15,10))

                qual,disc = plotall('{}/{}'.format(self.qa_dir,filestem))


            self.write_qa(filestem,qual,disc)
            if 'p' in qual:
                # I.E. if the data quality is good (or OK)
                print('Event rejected')
            else:
                self.accepted_i.append(i)
                print('Event accepted')

        #Slice the rows which correspond to good data and wirte them to a new pairs file
        print('{} accepted'.format(len(self.accepted_i)))
        accepted_pairs = self.pairs.iloc[self.accepted_i,:]
        #print(accepted_pairs)
        accepted_pairs.to_csv('{}/Accepted_SKS_SKKS.pairs'.format(self.result_path),sep=' ',index=False,mode='a+')

if __name__ == '__main__' :
    print('Hello World, this is quality_check.py. You are running me from the command line')

    pair_file = sys.argv[1]
    results_dir = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split'

    if len(sys.argv) is 1:
        print('No inputs detected. Please input them now')
        pair_file = input('Input pair filename \n >')

    Review = Inspecter(pair_file,results_dir)


    Review.run_inspection()