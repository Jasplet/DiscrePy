
#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.gridspec as gridspec
import obspy
from obspy import taup
import os
from subprocess import call
import shlex

def make_pairs(path,sdb_stem):
    """
    Function to take .sdb files for SKS and SKKS and generate the set of SKS-SKKS pais that exist in the datasetself.
    inputs
    - sdb_stem: [str] sdb file stem for the data that you want to process
    - path: [str] path to Sheba/Results directory that contains the sdb's that you want
    -

    """
    ## Set Path to Sheba Directory
    SKS_sdb = '{}_SKS_sheba_results.sdb'.format(sdb_stem)
    SKKS_sdb = '{}_SKKS_sheba_results.sdb'.format(sdb_stem)
    # First import the SKS and SKKS .sdb files (sdb = splitting database)
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    SKS = pandas.read_csv('{}/{}'.format(path,SKS_sdb),delim_whitespace=True,converters=date_time_convert)
    SKKS = pandas.read_csv('{}/{}'.format(path,SKKS_sdb),delim_whitespace=True,converters=date_time_convert)
    # Now the sdb files have been read as pandas dataframes, we can perform an inner join. This will return a single dataframe containing all rows from SKS and SKKS where
    # ['DATE','TIME','STAT','STLA','STLO','EVLA','EVLO','EVDP','DIST','AZI','BAZ'] are the same.
    SKS_SKKS_pair = pandas.merge(SKS,SKKS,on=['DATE','TIME','STAT','STLA','STLO','EVLA','EVLO','EVDP','DIST','AZI','BAZ'],how='inner')
    relabel = {'FAST_x':'FAST_SKS', 'DFAST_x': 'DFAST_SKS','TLAG_x':'TLAG_SKS','DTLAG_x':'DTLAG_SKS','SPOL_x':'SPOL_SKS','DSPOL_x':'DSPOL_SKS',
          'WBEG_x':'WBEG_SKS','WEND_x':'WEND_SKS','EIGORIG_x':'EIGORIG_SKS','EIGCORR_x':'EIGCORR_SKS','Q_x':'Q_SKS','SNR_x':'SNR_SKS','NDF_x':'NDF_SKS','FAST_y':'FAST_SKKS', 'DFAST_y': 'DFAST_SKKS',
          'TLAG_y':'TLAG_SKKS','DTLAG_y':'DTLAG_SKKS','SPOL_y':'SPOL_SKKS','DSPOL_y':'DSPOL_SKKS','WBEG_y':'WBEG_SKKS','WEND_y':'WEND_SKKS','EIGORIG_y':'EIGORIG_SKKS','EIGCORR_y':'EIGCORR_SKKS',
          'Q_y':'Q_SKKS','SNR_y':'SNR_SKKS','NDF_y':'NDF_SKKS'}
    # The dictionary relabels the other columns in the join so that we can more easily pick apart the SKS and SKKS results
    SKS_SKKS_pair.rename(relabel,axis='columns',inplace=True)
    # Sort the Pairs dataframe so the pairs are in chronological order (by origin time (DATE only))
    s = SKS_SKKS_pair.sort_values(by=['DATE'],ascending=True)
    # Save the dataframe to a new pairs file. The stem of the sdb files is used so that is it clear which .pairs file relates to which .sdb files
    s.to_csv('{}/{}_SKS_SKKS.pairs'.format(path,sdb_stem),sep=' ',index=False)


class Pairs:
    """
    Class to hold a .pairs Database (and Pierce Points) and generate a suite of useful plots based off the data.
    """
    def __init__(self,pf,Q_threshold=None,gcarc_threshold=None):

        self.pairs = pd.read_csv('{}'.format(pf),delim_whitespace=True,converters={'TIME': lambda x: str(x)})
        #Load raw data from provided .pairs file. This is going to be a hidden file as I will parse out useful columns to new attributes depending of provided kwargs

         #if kwargs are none:\

        if os.path.isfile('{}.pp'.format(pf[:-6])):
            #print(pf[:-6])
            self.pp = pd.read_csv('{}.pp'.format(pf[:-6]),delim_whitespace=True)
        else:
            print('Pierce Points file {}.pp doesnt not exist, calling pierce.sh'.format(pf[:-6]))
            p = call(shlex.split('/Users/ja17375/Shear_Wave_Splitting/Sheba/Programs/pierce.sh {}'.format(pf)))
            self.pp = self.pp = pd.read_csv('{}.pp'.format(pf[:-6]),delim_whitespace=True)
        # Load SDB and PP (pierce points data for a set of SKS-SKKS pairs)



    def match(self,path,filestem,sigma=2):
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

        outfile = open('{}/{}_matches.pairs'.format(path,filestem),'w+')
        outfile2 = open('{}/{}_diffs.pairs'.format(path,filestem),'w+')
        mspp1 = open('{}/{}_matches.mspp'.format(path,filestem),'w+')
        mspp2 = open('{}/{}_diffs.mspp'.format(path,filestem),'w+')
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

#####################################################################################################
# Top level where script is invoked from command line
if __name__ == '__main__':
    print('Hello World, I am SDB_analysis.py and now its time to do my thing')

    if len(sys.argv) < 3:
        print('Excuse me, but you appear to have forgotten to provide enough arguements, 2 arguements (path and sdb_stem) are needed. Lets do this now.')
        path = input('Input the path to the Sheba/Results directory that contains the sdb files \n > ')
        sdb_stem = input('Input the stem of the SKS and SKKS sdb files you want to analse \n > ')


    elif len(sys.argv == 3):
        path = sys.argv[1]
        sdb_stem = sys.argv[2]

    else:
        Warning('Too many arguements, script will not do anything else')

    if os.path.isfile('{}/{}'.format(path,sdb_stem)) is False:
        "Lets's find some pairs"
        make_pairs(path,sdb_stem)
    else:
        'Pairs already exist'

    # Now lets read the pairs
    pair_file = '{}/{}_SKS_SKKS.pairs'.format(path,sdb_stem)
    p = Pairs(pair_file)
    p.match(path,pair_file[:-6])

    print('End')
