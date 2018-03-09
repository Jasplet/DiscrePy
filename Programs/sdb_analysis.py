
#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import cartopy.crs as cart
#import cartopy
import matplotlib.gridspec as gridspec
import obspy
from obspy import taup



class SDB:
    """
    Class to hold a Splitting Database (and Pierce Points) and generate a suite of useful plots based off the data.
    """
    def __init__(self,sdb,Q_threshold=None,gcarc_threshold=None):

        self._raw = pd.read_csv('{}.sdb'.format(sdb),delim_whitespace=True)
        #Load raw data from provided sdb file. This is going to be a hidden file as I will parse out useful columns to new attributes depending of provided kwargs

         #if kwargs are none:
        self.sdb = self._raw
        self.pp = pd.read_csv('{}.pp'.format(sdb),delim_whitespace=True)
        ## Load SDB and PP (pierce points data for a set of SKS-SKKS pairs)

    def main(self):
        """
        Main Function of the Class
        """

    def match(self,sigma=2):
        """
        Funntion to see if the SKS and SKKS splititng measurements for a pair of measurements match within error

        Default error for this kind of anlysis is 2-sigma. Sheba returns 1 sigma so the DFAST and DTLAG need to be scaled appropriatly.
        """
        #First lets extract the raw values of the data that we need
        SKS_fast = self.data.FAST_SKS.values
        SKS_dfast = self.data.DFAST_SKS.values
        SKS_tlag = self.data.TLAG_SKS.values
        SKS_dtlag = self.data.DTLAG_SKS.values
        #Niw for SKKS
        SKKS_fast = self.data.FAST_SKKS.values
        SKKS_dfast = self.data.DFAST_SKKS.values
        SKKS_tlag = self.data.TLAG_SKKS.values
        SKKS_dtlag = self.data.DTLAG_SKKS.values
        # Now set the SKS and SKKS 2-sigma rnages
        lbf_SKS = SKS_fast - sigma*SKS_dfast
        ubf_SKS = SKS_fast + sigma*SKS_dfast
        lbf_SKKS = SKKS_fast - sigma*SKKS_dfast
        ubf_SKKS = SKKS_fast + sigma*SKKS_dfast

        lbt_SKS = SKS_tlag - sigma*SKS_dtlag
        ubt_SKS = SKS_tlag - sigma*SKS_dtlag
        lbt_SKKS = SKKS_tlag - sigma*SKKS_dtlag
        ubt_SKKS = SKKS_tlag - sigma*SKKS_dtlag
