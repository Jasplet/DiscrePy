
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
        SKS_fast = self.sdb.FAST_SKS.values
        SKS_dfast = self.sdb.DFAST_SKS.values
        SKS_tlag = self.sdb.TLAG_SKS.values
        SKS_dtlag = self.sdb.DTLAG_SKS.values
        #Niw for SKKS
        SKKS_fast = self.sdb.FAST_SKKS.values
        SKKS_dfast = self.sdb.DFAST_SKKS.values
        SKKS_tlag = self.sdb.TLAG_SKKS.values
        SKKS_dtlag = self.sdb.DTLAG_SKKS.values
        # Now set the SKS and SKKS 2-sigma rnages
        lbf_SKS = SKS_fast - sigma*SKS_dfast
        ubf_SKS = SKS_fast + sigma*SKS_dfast
        lbf_SKKS = SKKS_fast - sigma*SKKS_dfast
        ubf_SKKS = SKKS_fast + sigma*SKKS_dfast

        lbt_SKS = SKS_tlag - sigma*SKS_dtlag
        ubt_SKS = SKS_tlag - sigma*SKS_dtlag
        lbt_SKKS = SKKS_tlag - sigma*SKKS_dtlag
        ubt_SKKS = SKKS_tlag - sigma*SKKS_dtlag

        outfile = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/SKS_SKKS_matches.sdb','w+')
        outfile2 = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/SKS_SKKS_diffs.sdb','w+')
        outfile.write('DATE TIME STAT SKS_PP_LAT SKS_PP_LON SKKS_PP_LAT SKKS_PP_LON SKS_FAST SKS_DFAST SKS_TLAG SKS_DTLAG SKKS_FAST SKKS_DFAST SKKS_TLAG SKKS_DTLAG\n')
        outfile2.write('DATE TIME STAT SKS_PP_LAT SKS_PP_LON SKKS_PP_LAT SKKS_PP_LON SKS_FAST SKS_DFAST SKS_TLAG SKS_DTLAG SKKS_FAST SKKS_DFAST SKKS_TLAG SKKS_DTLAG\n')
        for i in enumerate(SKS_fast):
            date = self.sdb.DATE.values[i]
            time = self.sdb.TIME.values[i]
            stat = self.sdb.STAT.values[i]
            SKS_pp_lat = self.pp.lat_SKS.values[i]
            SKS_pp_lon = self.pp.lon_SKS.values[i]
            SKKS_pp_lat = self.pp.lat_SKKS.values[i]
            SKKS_pp_lon = self.pp.lon_SKKS.values[i]
            if (lbf_SKKS[i] =< SKS_fast[i] =< ubf_SKKS[i] ) and (lbt_SKKS[i] =< SKS_tlag[i] =< ubt_SKKS[i]):
                # Do the Fast and Tlag measured for SKS sit within the error bars for SKKS?
                if (lbf_SKS[i] <= SKKS_fast[i] =< ubf_SKS[i]) and (lbt_SKS[i] =< SKKS_tlag =< ubt_SKS[i]):
                    # Do the Fast and Tlag measured for SKKS sit within the 2-sigma error bars for SKS?

                    outfile.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
                else:
                    outfile2.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
            else:
                outfile2.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))

            #End of the matching If block

        outfile.close()
