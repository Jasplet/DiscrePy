#! /usr/bin/env python
####################################
#   Program: sac_download.py
####################################
#   Author: Joseph Asplet
####################################
#   This is a program to download seismograms for a given set
#   of events in the SAC format.
####################################
#   Standard Import Statements
####################################
import obspy as ob
import pandas as pd
import numpy as np
import obspy.core.utcdatetime
import os.path
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import (FDSNException)
from obspy.clients.fdsn.header import (FDSNNoDataException)
from obspy.core import AttribDict
from obspy.clients import iris
####################################
def main(event_list=None,batch=False):
    """
    Main Routine for this program. Holds highest level logic
    """
    if batch is True:
        # Reads our event list as a pandas datatframe
        df = pd.read_csv(eventlist,delim_whitespace=True)
        for station in df.STAT.unique():
#        for each unique station Code
            Instance = Downloader(df,station)


class Downloader:

    def __init__(self,df,station):

        self.station = station
        self.data = df[(df['STAT'] == station)]
        self.data = self.data.reset_index()
        del self.data['index']
#           Resets indexing of DataFrame


        try:
            print('Make /Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}'.format(station))
            os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}'.format(station))
        except FileExistsError:
            print('It already exists, Hooray! Less work for me!')
    #   Made
        self.outfile = open('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files/{}/{}_downloaded_streams.txt'.format(station,station),'w+')

        self.attempts = 0 #Counter for how many attempted downloads there were
        self.fdsnx = 0 #Counter for how many attempts hit a FDSNNoDataException
        self.dwn = 0 #Counter for how many events were downloaded
        self.ex = 0 #Counter for how many event already exist in filesystem and therefore werent downloaded
        self.ts = 0 #Counter for events who;s traces are too short.


def download():
    """
    Function to hold download requests
    """

# PsuedoCode - For developoment

# read list of events as daatframe
#    for unique stations:
#       get network/location code information
#       for events in df.stat=stat:
#           get accurate event time
#               for BHN,BHE.BHZ:
#                   request download of trace of set length
#                   save to holder
#                   read holder and add desired sac headers
#
