
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
import shlex



class Pairs:
    """
    Class to hold a .pairs Database (and Pierce Points) and generate a suite of useful plots based off the data.
    """
    def __init__(self,pf,Q_threshold=None,gcarc_threshold=None):

        self.pairs = pd.read_csv('{}'.format(pf),delim_whitespace=True,converters={'TIME': lambda x: str(x)})
        #Load raw data from provided sdb file. This is going to be a hidden file as I will parse out useful columns to new attributes depending of provided kwargs

         #if kwargs are none:\

        if os.path.isfile('{}.pp'.format(pf[:-6])):
            #print(pf[:-6])
            self.pp = pd.read_csv('{}.pp'.format(pf[:-6]),delim_whitespace=True)
        else:
            print('Pierce Points file {}.pp doesnt not exist, calling pierce.sh'.format(pf[:-6]))
            p = call(shlex.split('/Users/ja17375/Shear_Wave_Splitting/Sheba/Programs/pierce.sh {}'.format(pf)))
            self.pp = self.pp = pd.read_csv('{}.pp'.format(pf[:-6]),delim_whitespace=True)
        # Load SDB and PP (pierce points data for a set of SKS-SKKS pairs)

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
        SKS_fast = self.pairs.FAST_SKS.values
        SKS_dfast = self.pairs.DFAST_SKS.values
        SKS_tlag = self.pairs.TLAG_SKS.values
        SKS_dtlag = self.pairs.DTLAG_SKS.values
        #Niw for SKKS
        SKKS_fast = self.pairs.FAST_SKKS.values
        SKKS_dfast = self.pairs.DFAST_SKKS.values
        SKKS_tlag = self.pairs.TLAG_SKKS.values
        SKKS_dtlag = self.pairs.DTLAG_SKKS.values
        # Now set the SKS and SKKS 2-sigma rnages
        lbf_SKS = SKS_fast - sigma*SKS_dfast
        ubf_SKS = SKS_fast + sigma*SKS_dfast
        lbf_SKKS = SKKS_fast - sigma*SKKS_dfast
        ubf_SKKS = SKKS_fast + sigma*SKKS_dfast

        lbt_SKS = SKS_tlag - sigma*SKS_dtlag
        ubt_SKS = SKS_tlag + sigma*SKS_dtlag
        lbt_SKKS = SKKS_tlag - sigma*SKKS_dtlag
        ubt_SKKS = SKKS_tlag + sigma*SKKS_dtlag

        outfile = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}_matches.pairs'.format(file),'w+')
        outfile2 = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}_diffs.pairs'.format(file),'w+')
        mspp1 = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}_matches.mspp'.format(file),'w+')
        mspp2 = open('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}_diffs.mspp'.format(file),'w+')
        outfile.write('DATE TIME STAT STLA STLO EVLA EVLO SKS_PP_LAT SKS_PP_LON SKKS_PP_LAT SKKS_PP_LON SKS_FAST SKS_DFAST SKS_TLAG SKS_DTLAG SKKS_FAST SKKS_DFAST SKKS_TLAG SKKS_DTLAG\n')
        outfile2.write('DATE TIME STAT STLA STLO EVLA EVLO SKS_PP_LAT SKS_PP_LON SKKS_PP_LAT SKKS_PP_LON SKS_FAST SKS_DFAST SKS_TLAG SKS_DTLAG SKKS_FAST SKKS_DFAST SKKS_TLAG SKKS_DTLAG\n')
        for i,value in enumerate(SKS_fast):
            date = self.pairs.DATE.values[i]
            time = self.pairs.TIME.values[i]
            stat = self.pairs.STAT.values[i]
            evla = self.pairs.EVLA.values[i]
            evlo = self.pairs.EVLO.values[i]
            stla = self.pairs.STLA.values[i]
            stlo = self.pairs.STLO.values[i]
            SKS_pp_lat = self.pp.lat_SKS.values[i]
            SKS_pp_lon = self.pp.lon_SKS.values[i]
            SKKS_pp_lat = self.pp.lat_SKKS.values[i]
            SKKS_pp_lon = self.pp.lon_SKKS.values[i]
            print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            #print('{} <= {} <= {} and {} <= {} <= {}'.format(lbf_SKKS[i],SKS_fast[i],ubf_SKKS[i],lbt_SKKS[i],SKS_tlag[i],ubt_SKKS[i]))
            if (lbf_SKKS[i] <= ubf_SKS[i]) or (lbf_SKS[i] <= ubf_SKKS[i]):
                # Do the Fast and Tlag measured for SKS sit within the error bars for SKKS?
                #print('{} <= {} <= {} and {} <= {} <= {}'.format(lbf_SKKS[i],SKS_fast[i],ubf_SKKS[i],lbt_SKKS[i],SKS_tlag[i],ubt_SKKS[i]))
                if (lbt_SKKS[i] <= ubt_SKS[i]) and (lbt_SKS[i] <= ubt_SKKS[i]):
                    # Do the Fast and Tlag measured for SKKS sit within the 2-sigma error bars for SKS?
                    #print('{} <= {} <= {} and {} <= {} <= {}'.format(lbf_SKS[i],SKKS_fast[i],ubf_SKS[i],lbt_SKS[i],SKKS_tlag[i],ubt_SKS[i]))
                    outfile.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,stla,stlo,evla,evlo,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
                    mspp1.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))
                else:
                    outfile2.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,stla,stlo,evla,evlo,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
                    mspp2.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))
            else:
                outfile2.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(date,time,stat,stla,stlo,evla,evlo,SKS_pp_lat,SKS_pp_lon,SKKS_pp_lat,SKKS_pp_lon,SKS_fast[i],SKS_dfast[i],SKS_tlag[i],SKS_dtlag[i],SKKS_fast[i],SKKS_dfast[i],SKKS_tlag[i],SKKS_dtlag[i]))
                mspp2.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))
            #End of the matching If block

        outfile.close()
        outfile2.close()
        mspp1.close()
        mspp2.close()

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
        sacfiles = list(map(mk_sacfiles_mypath,self.pairs.STAT,self.pairs.DATE,self.pairs.TIME))

        for file in sacfiles:
            cp_sacfile(file,mypath,outdir)
