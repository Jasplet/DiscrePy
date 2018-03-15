
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
import os
from subprocess import call
from multiprocessing import Pool



class SDB:
    """
    Class to hold a Splitting Database (and Pierce Points) and generate a suite of useful plots based off the data.
    """
    def __init__(self,sdb,Q_threshold=None,gcarc_threshold=None):

        self._raw = pd.read_csv('{}.sdb'.format(sdb),delim_whitespace=True,converters={'TIME': lambda x: str(x)})
        #Load raw data from provided sdb file. This is going to be a hidden file as I will parse out useful columns to new attributes depending of provided kwargs

         #if kwargs are none:
        self.sdb = self._raw
        self.pp = pd.read_csv('{}.pp'.format(sdb),delim_whitespace=True)
        ## Load SDB and PP (pierce points data for a set of SKS-SKKS pairs)

    def main(self):
        """
        Main Function of the Class
        """

    def match(self,file,sigma=2):
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
        ubt_SKS = SKS_tlag + sigma*SKS_dtlag
        lbt_SKKS = SKKS_tlag - sigma*SKKS_dtlag
        ubt_SKKS = SKKS_tlag + sigma*SKKS_dtlag

        outfile = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/{}_matches.sdb'.format(file),'w+')
        outfile2 = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/{}_diffs.sdb'.format(file),'w+')
        outfile.write('DATE TIME STAT SKS_PP_LAT SKS_PP_LON SKKS_PP_LAT SKKS_PP_LON SKS_FAST SKS_DFAST SKS_TLAG SKS_DTLAG SKKS_FAST SKKS_DFAST SKKS_TLAG SKKS_DTLAG\n')
        outfile2.write('DATE TIME STAT SKS_PP_LAT SKS_PP_LON SKKS_PP_LAT SKKS_PP_LON SKS_FAST SKS_DFAST SKS_TLAG SKS_DTLAG SKKS_FAST SKKS_DFAST SKKS_TLAG SKKS_DTLAG\n')
        for i,value in enumerate(SKS_fast):
            date = self.sdb.DATE.values[i]
            time = self.sdb.TIME.values[i]
            stat = self.sdb.STAT.values[i]
            SKS_pp_lat = self.pp.lat_SKS.values[i]
            SKS_pp_lon = self.pp.lon_SKS.values[i]
            SKKS_pp_lat = self.pp.lat_SKKS.values[i]
            SKKS_pp_lon = self.pp.lon_SKKS.values[i]
            #print('{} <= {} <= {} and {} <= {} <= {}'.format(lbf_SKKS[i],SKS_fast[i],ubf_SKKS[i],lbt_SKKS[i],SKS_tlag[i],ubt_SKKS[i]))
            if (lbf_SKKS[i] <= SKS_fast[i] <= ubf_SKKS[i] ) and (lbt_SKKS[i] <= SKS_tlag[i] <= ubt_SKKS[i]):
                # Do the Fast and Tlag measured for SKS sit within the error bars for SKKS?
                #print('{} <= {} <= {} and {} <= {} <= {}'.format(lbf_SKKS[i],SKS_fast[i],ubf_SKKS[i],lbt_SKKS[i],SKS_tlag[i],ubt_SKKS[i]))
                if (lbf_SKS[i] <= SKKS_fast[i] <= ubf_SKS[i]) and (lbt_SKS[i] <= SKKS_tlag[i] <= ubt_SKS[i]):
                    # Do the Fast and Tlag measured for SKKS sit within the 2-sigma error bars for SKS?
                    print('{} <= {} <= {} and {} <= {} <= {}'.format(lbf_SKS[i],SKKS_fast[i],ubf_SKS[i],lbt_SKS[i],SKKS_tlag[i],ubt_SKS[i]))
                    outfile.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
                else:
                    outfile2.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
            else:
                outfile2.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))

            #End of the matching If block

        outfile.close()



    def package(self,outdir,mypath='/Users/ja17375/Shear_Wave_Splitting/Sheba/SAC'):
        '''
        Function to package the SAC files for the observations made into a new directory (for sharing data and results)
        '''
        if os.path.isdir('{}{}'.format(mypath,outdir)) is False:
            #If outdir doesnt exist then Make it ad underlying structure (an additional SAC directory )
            os.makedirs('{}{}/SAC'.format(mypath,outdir))

        def mk_sacfiles(path,stat,date,time):
            '''Funciton to generate sacfile names'''
            #path='/Users/ja17375/Shear_Wave_Splitting/Sheba/SAC'
            return '{}/{}/{}_{}_{}*BH?.sac'.format(path,stat,stat,date,time)

        def cp_sacfile(sacfile,path,outdir):
            '''Function to copy sacfile to the output directory
            Path - path to the output directorys
            Outdir - output directory name
            '''
            call('cp {} {}/{}/SAC/'.format(sacfile,path,outdir),shell=True)
            print('cp {} {}/{}/SAC/'.format(sacfile,path,outdir))

        mk_sacfiles_mypath = lambda stat,date,time : mk_sacfiles('/Users/ja17375/Shear_Wave_Splitting/Data/SAC_files',stat,date,time)
        cp_sacfiles_mypath = lambda sacfile: cp_sacfile(sacfile,mypath,outdir)
        sacfiles = list(map(mk_sacfiles_mypath,self.sdb.STAT,self.sdb.DATE,self.sdb.TIME))

        for file in sacfiles:
            cp_sacfile(file,mypath,outdir)
