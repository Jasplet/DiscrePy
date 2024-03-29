#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import sys
import os
import shlex
from subprocess import call
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib
import numpy as np
from stack import Stacker
from time import ctime
from glob import glob
from l2stats import ftest
import obspy
import obspy.taup


class Builder:
    """
    Class to construct a .pairs Database (and Pierce Points), calc D_SI values and stack L2 surfaces.

    - sdb_stem: [str] sdb file stem for the data that you want to process
    - path: [str] path to Sheba/Results directory that contains the sdb's that you want
    - RunDir: [str] path to the Run Directory than Contains the .lamR surfaces to stack
    """
    def __init__(self,p,RunDir,sdb_stem,snr=5.0,syn=False):

        self.path = p # Path to Results directory
        self.path_stk = '/Users/ja17375/DiscrePy/Sheba/Runs/{}'.format(RunDir)
        self.sdb_stem = sdb_stem
        self.fpath =     '{}/{}_{:02d}.pairs'.format(self.path,self.sdb_stem,int(snr))
        self.lam2_bar = [ ] # Variabel to hold lambda2 values returned from stacking routine
        self.lam2alpha_sks = [ ] # Variable to hold list of lam2alpha values for SKS
        self.lam2alpha_skks = [ ] # Variable to holf lam2alpha values for SKKS
        self.lam2_sum = [ ]
        self.snr = snr # Holds SNR cutoff (defaul is 5)
        self.syn=syn
        #if kwargs are none:\

    def run(self):
        '''Function to build the pairs file'''
        # First match the SKS and SKKS
        start = ctime()
        print('Making Pairs')
        self.make_pairs()
        # Apply a quick Signal to Noise test to get rid of the rreally bad data
        self.P = self.QA_tests() # Overwrite self.P weith only accepted events
        self.write_out(self.P,name='{}_{:02d}.pairs'.format(self.sdb_stem,int(self.snr)))
        # Write initial pairs file so we can make piercepoints
        # Next generate the piercepoints and add them to the df
        print('Adding PiercePoints')
        self.gen_pp()
        self.add_pp()
        # Now calculate the d_SI vlaues
        print('Calculate difference in splitting intensity')
        self.add_DSI()
        # Finally stack the lamR surfaces to determine the lam2 values
        if self.syn is False:
            ''' I.e Are we using real data? Because of the set up of the synthetics (and I am only using 6 pairs) I will stack these more manually and I dont want to reject any '''
            print('Calculate lambda 2 values')
            self.add_lam2()
            #Now test for matching and disrecpent pairs
            print('Apply 2-sigma test for discrepancy')
            # self.match_sig2()
            print(r'Apply $\bar{\lambda_2}$ & $\Delta SI$ test')
            self.match_l2()
            # print('{} pairs'.format(len(self.P)))

        # And save the result
        self.write_out(self.P,'{}_{:02d}.pairs'.format(self.sdb_stem,int(self.snr)))
        end = ctime()
        print('END. start {}, end {}'.format(start,end))



    def gen_pp(self):
        ''' Fucntion to test for whether the .pp file exists and if not call TauP to generate it and the corresponding 
        files '''

        print('Looking for pp file {}.pp'.format(self.fpath.split('.')[0]))
        if os.path.isfile('{}.pp'.format(self.fpath.split('.')[0])):
            #print(pf[:-6])
            print('Exists')
            self.pp = pd.read_csv('{}.pp'.format(self.fpath.split('.')[0]),delim_whitespace=True)
        else:
            print('Pierce Points file {}.pp doesnt not exist, calling pierce.sh'.format(self.fpath.split('.')[0]))
            p = call(shlex.split('/Users/ja17375/DiscrePy/Sheba/Programs/pierce.sh {}'.format(self.fpath)))
            self.pp = pd.read_csv('{}.pp'.format(self.fpath.split('.')[0]),delim_whitespace=True)
        # Load SDB and PP (pierce points data for a set of SKS-SKKS pairs)

        if os.path.isfile('{}.mspp'.format(self.fpath.split('.')[0])) is False:
            print('{}.mspp does not exist, creating'.format(self.fpath.split('.')[0]))
            with open('{}.mspp'.format(self.fpath.split('.')[0]),'w+') as writer:
                for i,row in enumerate(self.pp.index):
                    writer.write('> \n {} {} \n {} {} \n'.format(self.pp.lon_SKS[i],self.pp.lat_SKS[i],self.pp.lon_SKKS[i],self.pp.lat_SKKS[i]))

    def make_pairs(self):
        """
        Function to take .sdb files for SKS and SKKS and generate the set of SKS-SKKS pais that exist in the dataset.
        inputs

        """
        ## Set Path to Sheba Directory
        SKS_sdb = '{}_SKS_sheba_results.sdb'.format(self.sdb_stem)
        SKKS_sdb = '{}_SKKS_sheba_results.sdb'.format(self.sdb_stem)
        # For Synthetics I didnt bother putting in the sheba results part so use these
        # SKS_sdb = '{}_SKS.sdb'.format(self.sdb_stem)
        # SKKS_sdb = '{}_SKKS.sdb'.format(self.sdb_stem)
        # First import the SKS and SKKS .sdb files (sdb = splitting database)
        date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
        SKS = pd.read_csv('{}/{}'.format(self.path,SKS_sdb),delim_whitespace=True,converters=date_time_convert)
        SKKS = pd.read_csv('{}/{}'.format(self.path,SKKS_sdb),delim_whitespace=True,converters=date_time_convert)
        # Now the sdb files have been read as pd dataframes, we can perform an inner join. This will return a single dataframe containing all rows from SKS and SKKS where
        # ['DATE','TIME','STAT','STLA','STLO','EVLA','EVLO','EVDP','DIST','AZI','BAZ'] are the same.
        SKS_SKKS_pair = pd.merge(SKS,SKKS,on=['DATE','TIME','STAT','STLA','STLO','EVLA','EVLO','EVDP','DIST','AZI','BAZ'],how='inner')
        relabel = {'FAST_x':'FAST_SKS', 'DFAST_x': 'DFAST_SKS','TLAG_x':'TLAG_SKS','DTLAG_x':'DTLAG_SKS','SPOL_x':'SPOL_SKS','DSPOL_x':'DSPOL_SKS',
              'WBEG_x':'WBEG_SKS','WEND_x':'WEND_SKS','EIGORIG_x':'EIGORIG_SKS','EIGCORR_x':'EIGCORR_SKS','Q_x':'Q_SKS','SNR_x':'SNR_SKS','NDF_x':'NDF_SKS',
              'SI(Pr)_x':'SI(Pr)_SKS', 'SI(Pa)_x':'SI(Pa)_SKS','FAST_y':'FAST_SKKS', 'DFAST_y': 'DFAST_SKKS','TLAG_y':'TLAG_SKKS','DTLAG_y':'DTLAG_SKKS',
              'SPOL_y':'SPOL_SKKS','DSPOL_y':'DSPOL_SKKS','WBEG_y':'WBEG_SKKS','WEND_y':'WEND_SKKS','EIGORIG_y':'EIGORIG_SKKS','EIGCORR_y':'EIGCORR_SKKS',
              'Q_y':'Q_SKKS','SNR_y':'SNR_SKKS','NDF_y':'NDF_SKKS','SI(Pr)_y':'SI(Pr)_SKKS', 'SI(Pa)_y':'SI(Pa)_SKKS'}
        # The dictionary relabels the other columns in the join so that we can more easily pick apart the SKS and SKKS results
        SKS_SKKS_pair.rename(relabel,axis='columns',inplace=True)
        # Sort the Pairs dataframe so the pairs are in chronological order (by origin time (DATE only))
        self.P = SKS_SKKS_pair.sort_values(by=['DATE'],ascending=True)

    def write_out(self,df,name):
        print('Writing to {}'.format(name))
        df.to_csv('{}/{}'.format(self.path,name),sep=' ',index=False)

    def add_DSI(self):
        '''Calculate the difference in Splitting Intensity for each pair and add it to dataframe'''
        si_pr_sks = self.P['SI(Pr)_SKS']
        si_pr_skks = self.P['SI(Pr)_SKKS']
        si_pa_sks = self.P['SI(Pa)_SKS']
        si_pa_skks = self.P['SI(Pa)_SKKS']
        d_si_pr = np.abs(si_pr_sks-si_pr_skks)
        d_si_pa = np.abs(si_pa_sks-si_pa_skks)
        # d = {'D_SI_Pr': d_si_pr,'D_SI_Pa':d_si_pa}
        d = {'SI_Pr_sks': si_pr_sks, 'SI_Pr_skks': si_pr_skks,'SI_Pa_sks': si_pa_sks,'SI_Pa_skks': si_pa_skks,
             'D_SI_Pr': d_si_pr,'D_SI_Pa':d_si_pa} # Special version which also adds raw splitting intensity
        ddf = pd.DataFrame(d)
        self.P[['SI_Pr_sks','SI_Pr_skks','SI_Pa_sks','SI_Pa_skks','D_SI_Pr','D_SI_Pa']] = ddf
        #Delete SI cols as we dont need them any more ?
        # del self.P['INTENS_x']
        # del self.P['INTENS_y']

    def pair_stack(self):
        ''' Runs Stacker for all the desired pairs (a .pairs file)'''

        ext ='lamR'
        print('Stacking')
        rd = self.path_stk.split('/')[-1]
        out = '/Users/ja17375/DiscrePy/Sheba/Results/Combined/{}/Stacks'.format(rd) # For Filt 03/05 casesed need to hardcode in Combined/ directory
        if os.path.isdir(out) is False:
            print('{} does not exist, creating'.format(out))
            os.mkdir('/Users/ja17375/DiscrePy/Sheba/Results/Combined/{}/Stacks'.format(rd))

        for i,f in enumerate(self.P.DATE.values):
            # print(len(self.P))
            # print('It {}, time is {} '.format(i,str(datetime.now())))
            # First get the right DATE,TIME and STATION
            date,time,stat = self.P.DATE[i], self.P.TIME[i], self.P.STAT[i]
            fstem = '{}_{}_{}'.format(stat,date,time)

            lam2_stem = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path_stk,stat,fstem,ext))
            # print(lam2_stem)
            print('{}/{}/SKS/{}??_SKS.lamR'.format(self.path_stk,stat,fstem))
            print('{}/{}/SKKS/{}??_SKKS.{}'.format(self.path_stk,stat,fstem,ext))
            if len(lam2_stem) != 0:
                # I.e if glob has managed to find the sks lam2 surface file
                sks_lam2 = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path_stk,stat,fstem,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.{}'.format(self.path_stk,stat,fstem,ext))[0]
                Stk = Stacker(sks_lam2,skks_lam2,out)
                self.lam2_bar.append(Stk.lam2_bar)
                ndf_sks, ndf_skks = self.P.NDF_SKS[i],self.P.NDF_SKKS[i]
                self.lam2alpha_sks.append(ftest(Stk.lam2_sks,ndf_sks)) # Calc Lam2 Alpha and append it to list for SKS and SKKS
                self.lam2alpha_skks.append(ftest(Stk.lam2_skks,ndf_skks))
            elif len(glob('{}/{}/SKS/{}_{}_*_SKS.{}'.format(self.path_stk,stat,stat,date,ext))) !=0 :
                fstem2 = '{}_{}'.format(stat,date)
                # print('fstem2')
                # print('{}/{}/SKS/{}_*_SKS.{}'.format(self.path_stk,stat,fstem2,ext))
                sks_lam2 = glob('{}/{}/SKS/{}_*_SKS.{}'.format(self.path_stk,stat,fstem2,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}_*_SKKS.{}'.format(self.path_stk,stat,fstem2,ext))[0]
                # Now for a sanity check
                if (len(sks_lam2) != 0) or (len(skks_lam2) != 0):
                    Stk = Stacker(sks_lam2,skks_lam2,out)
                    self.lam2_bar.append(Stk.lam2_bar)
                    ndf_sks, ndf_skks = self.P.NDF_SKS[i],self.P.NDF_SKKS[i]
                    self.lam2alpha_sks.append(ftest(Stk.lam2_sks,ndf_sks)) # Calc Lam2 Alpha and append it to list for SKS and SKKS
                    self.lam2alpha_skks.append(ftest(Stk.lam2_skks,ndf_skks))
                else:
                    print('lam2 surfaces cannot be found, skipping')
                    pass
            else:
                print('Glob returns no files, see:')
                print(glob('{}/{}/SKS/{}_{}_*_SKS.{}'))
            # Now calulate LAM2ALPHA  for SKS and SKKS


    def add_lam2(self):
        '''
        Stack the associated raw lambda 2 surfaces (as output by sheba) for each SKS SKKS pair and find min value
        '''

        self.pair_stack()
        lam2_sum = [self.lam2alpha_sks[i] + self.lam2alpha_skks[i] for i in range(0,len(self.lam2alpha_sks))]
        print(len(self.lam2_bar),len(lam2_sum),len(self.lam2alpha_sks),len(self.lam2alpha_skks))
        l2df = {'LAM2_BAR' : self.lam2_bar,  'LAM2_SUM' : lam2_sum, 'LAM2A_SKS' : self.lam2alpha_sks, 'LAM2A_SKKS' : self.lam2alpha_skks}
        ldf = pd.DataFrame(l2df)
        self.P[['LAM2_BAR','LAM2_SUM','LAM2A_SKS','LAM2A_SKKS']] = ldf

    def add_pp(self):
        '''Adds piercepoints to .pairs file'''
        if len(self.pp) == len(self.P):
            print(len(self.pp) ,len(self.P))
            self.P['SKS_PP_LAT'] = self.pp.lat_SKS
            self.P['SKS_PP_LON'] = self.pp.lon_SKS
            self.P['SKKS_PP_LAT'] = self.pp.lat_SKKS
            self.P['SKKS_PP_LON'] = self.pp.lon_SKKS

            print('Added piercepoints')
        else:
            print('Dimension mismatch, not adding piercepoints')
    def match_sig2(self,sigma=2):
        """
        Funntion to see if the SKS and SKKS splititng measurements for a pair of measurements match within error
        file - the filename you want for the output files
        Default error for this kind of anlysis is 2-sigma. Sheba returns 1 sigma so the DFAST and DTLAG need to be scaled appropriatly.
        """
        fstem = self.path.split('.')[0]# Split('.')[0] takes off the extension
        # First Add piercepoints to pairas files
        # self.add_pp()
        # Set the SKS and SKKS 2-sigma rnages
        lbf_SKS = self.P.FAST_SKS - sigma*self.P.DFAST_SKS
        ubf_SKS = self.P.FAST_SKS + sigma*self.P.DFAST_SKS
        lbf_SKKS = self.P.FAST_SKKS - sigma*self.P.DFAST_SKKS
        ubf_SKKS = self.P.FAST_SKKS + sigma*self.P.DFAST_SKKS

        lbt_SKS = self.P.TLAG_SKS - sigma*self.P.DTLAG_SKS
        ubt_SKS = self.P.TLAG_SKS + sigma*self.P.DTLAG_SKS
        lbt_SKKS = self.P.TLAG_SKKS - sigma*self.P.DTLAG_SKKS
        ubt_SKKS = self.P.TLAG_SKKS+ sigma*self.P.DTLAG_SKKS
        # Now set the logic for the tests
        fast_test = (lbf_SKKS <= ubf_SKS) & (lbf_SKS <= ubf_SKKS)
        lag_test = (lbt_SKKS <= ubt_SKS) & (lbt_SKS <= ubt_SKKS)
        # Q_test = (self.P.Q_SKS <= -0.5 | self.P.Q_SKKS <= -0.5) & (self.P.Q_SKS >= 0.5 | self.P.Q_SKKS >= 0.5)

        match  = self.P[(fast_test & lag_test)] # Test for pairs that match within the given sigma range
        diff =  self.P.drop(index=match.index) # Remove matching pairs from original df to get the different pairs.
        # Write out matching and discepent dataframes
        match.to_csv('{}/{}_{:02d}_matches.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        diff.to_csv('{}/{}_{:02d}_diffs.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        # Open up mspp files
        print('Writing to {}/{}_{:02d}'.format(self.path,self.sdb_stem,int(self.snr)))
        mspp_match = open('{}/{}_{:02d}_matches.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_diff = open('{}/{}_{:02d}_diffs.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')

        for i,index in enumerate(diff.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_diff.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        for i,index in enumerate(match.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_match.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        mspp_diff.close()
        mspp_match.close()

    def match_l2(self):
        """
        Function to test if a SK(K)S pair are matching or discrepant. This test is based on lb2 and dSI (see synthetics_stakcs notes)
        This function assumes that lam2 and dSI have already been calculated and added to the pairs file. It also assumes that you
        want to keep the same file names but just add _match_l2 and _diff_l2self.
        The lam2 threshold used is 0.01 and the dSI threshold is 0.15
        """

        fstem = self.path.split('.')[0]# Split('.')[0] takes off the extension
        P = self.P.copy()
        # print(P.Q_SKS)
        # print(P.Q_SKKS)
        # Now apply the test to find the discrepant pairs, by definition the remainder must by the matches
        # First find pairs that are split or not
        uID = P[((P.Q_SKS > -0.7) & (P.Q_SKS < 0.5)) | ((P.Q_SKKS > -0.7) & (P.Q_SKKS < 0.5))]
        P.drop(uID.index)
        null_pairs  = P[((P.Q_SKS <= -0.7) & (P.Q_SKKS <= -0.7))] # Pairs where both phases are nulls (according to Q), auto classify as matching
        null_split_pair = P[(((P.Q_SKS <= -0.7) & (P.Q_SKKS >= 0.5)) | ((P.Q_SKS >= 0.5) & (P.Q_SKKS <= -0.7)))] # Test for pairs with 1 null 1 split, discrepant by definition
        splits = P[((P.Q_SKS > 0.5) & (P.Q_SKKS > 0.5 ))] # Test for pairs whjere both phases are split
        t_l2_splits = splits.LAM2_SUM # Lam2_SUm is now the sum of the two 95% confidence levels
        t_l2_ns = null_split_pair.LAM2_SUM
        t_dSI = 0.4 # Threshold of 0.4 taken from Deng et al (2017)
        diff= splits[(splits.LAM2_BAR > t_l2_splits) & (splits.D_SI_Pr > t_dSI)] #| (splits.D_SI > t_dSI))] # Apply tests for discrepant splitting
        match = splits[(splits.LAM2_BAR <= t_l2_splits) | (splits.D_SI_Pr <= t_dSI)] # If the pair fails either lam2 or dSI test then we call it matching
        diff_dsi = splits[(splits.D_SI_Pr > 0.4)]
        match_dsi = splits[(splits.D_SI_Pr <= 0.4)]
        ns_diff = null_split_pair[(null_split_pair.LAM2_BAR > t_l2_ns) & (null_split_pair.D_SI_Pr > t_dSI)]
        ns_match = null_split_pair[(null_split_pair.LAM2_BAR <= t_l2_ns) | (null_split_pair.D_SI_Pr <= t_dSI)]

        print(len(self.P))
        print('There are {} split pairs. {} are matches and {} are discrepant!'.format(len(splits),len(match),len(diff)))
        print('There are {} null-split pairs. {} are matches and {} are discrepant!'.format(len(null_split_pair),len(ns_match),len(ns_diff)))

        test = len(uID) + len(null_pairs) + len(null_split_pair) + len(diff) + len(match)
        print(test)
        # Now combined matching and discrepant pairs together
        null_pairs.to_csv('{}/{}_{:02d}_nulls.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        ns_match.to_csv('{}/{}_{:02d}_matches_null_split.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        ns_diff.to_csv('{}/{}_{:02d}_diffs_null_split.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        match.to_csv('{}/{}_{:02d}_matches.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        diff.to_csv('{}/{}_{:02d}_diffs.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        match_dsi.to_csv('{}/{}_{:02d}_matches_dsi.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        diff_dsi.to_csv('{}/{}_{:02d}_diffs_dsi.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        uID.to_csv('{}/{}_{:02d}_uID_l2.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        # Open up mspp files
        print('Writing to {}/{}_{:02d}'.format(self.path,self.sdb_stem,int(self.snr)))
        mspp_match = open('{}/{}_{:02d}_matches.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_diff = open('{}/{}_{:02d}_diffs.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        # mspp_match_dsi = open('{}/{}_{:02d}_matches_dsi.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        # mspp_diff_dsi = open('{}/{}_{:02d}_diffs_dsi.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_uID = open('{}/{}_{:02d}_uID_l2.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_null_pairs = open('{}/{}_{:02d}_nulls.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_null_split_match = open('{}/{}_{:02d}_matches_null_split.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_null_split_diff = open('{}/{}_{:02d}_diffs_null_split.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')


        for i,index in enumerate(diff.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_diff.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        for i,index in enumerate(null_pairs.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_null_pairs.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        for i,index in enumerate( ns_match.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_null_split_match.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        for i,index in enumerate( ns_diff.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_null_split_diff.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        for i,index in enumerate(match.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_match.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        for i,index in enumerate(uID.index):
            SKS_pp_lat = self.pp.lat_SKS.values[index]
            SKS_pp_lon = self.pp.lon_SKS.values[index]
            SKKS_pp_lat = self.pp.lat_SKKS.values[index]
            SKKS_pp_lon = self.pp.lon_SKKS.values[index]
            #print(i,date,stat,evla,evlo,stla,stlo,SKS_pp_lat,SKS_pp_lon)
            mspp_uID.write('> \n {} {} \n {} {} \n'.format(SKS_pp_lon,SKS_pp_lat,SKKS_pp_lon,SKKS_pp_lat))

        mspp_uID.close()
        mspp_diff.close()
        mspp_match.close()
        mspp_null_pairs.close()
        mspp_null_split_match.close()
        mspp_null_split_diff.close()

    def QA_tests(self):
        ''' Run some basic QA tests on SNR and SPOL/BAZ. We test to see if signal to noise is greater than 5 (which our synthetics work has shown to be the point where one can
        reasonable expect sheba to converge toward the true solution for phi,dt) and SPOL/BAZ to avoid beamformin imssues etc. SPOL/BAZ failuers are output for potential future investigation'''
        t = self.snr
        print(t)
        self.accepted_i = [ ]
        self.snr_fail = [ ]
        self.baz_spol_fail = [ ]
        print('There are {} pairs pre-SNR < 2 test'.format(len(self.P)))
        for i,row in self.P.iterrows():
            if row.SNR_SKS > self.snr and row.SNR_SKKS > self.snr:
                #Test to see if Signal-to-Noise is too high
                if abs(row.BAZ%180 - row.SPOL_SKS%180) <  10 and abs(row.BAZ%180 - row.SPOL_SKKS%180) < 10:
                # Test to see if there is a asignficiant difference between source polarisation and back azimuth abs
                    print('Event accepted')
                    self.accepted_i.append(i)
                else:
                    print('SPOL-BAZ difference greater than 10, reject')
                    self.baz_spol_fail.append(i)

            else:
                print('SNR for SKS or SKKS less than {:02}, auto-reject'.format(self.snr))
                self.snr_fail.append(i)


        print('{} accepted, {} rejected for SNR, {} rejected for BAZ-SPOL'.format(len(self.accepted_i),len(self.snr_fail),len(self.baz_spol_fail)))#,len(self.baz_spol_fail)))
        drop = self.snr_fail #+ self.baz_spol_fail
        accepted_pairs_snr = self.P.drop(drop)
        accepted_pairs = accepted_pairs_snr.drop(self.baz_spol_fail)
        rejected_baz_spol = accepted_pairs_snr.drop(self.accepted_i)
        rejected_baz_spol.to_csv('{}/{}_{:02d}_rejected_baz_spol.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        accepted_pairs.index = pd.RangeIndex(len(accepted_pairs.index))
        return accepted_pairs


class Pairs:
    """ A class to actually hold a Pairs (read as a pandas dataframe) and then do useful stuff with it """

    def __init__(self,df=None,file=False,fname=None,syn=False,syn1=None,syn2=None,synstk=None):
        '''
        df - a pandas dataframe contaiing the pairs (implicitly assumed that this df has been made using builder)
        file [bool] - T/F flag for if you want to read in a already existing pairs file. This option is easier if you have not freshly created it as the correct converters to preserve
                    leading zeros in date and time will be applied
        fname [str] - Full path to (and including) the pairs file your want to read in.
        '''
        self.rpath = '/Users/ja17375/DiscrePy/'
        self.syn = syn # [Bool] - is data synthetics or not?
        if self.syn == True:
            self.syn1 = syn1
            self.syn2 = syn2
            self.synstk = synstk

        if file is True:
            print('Read file option specified')
            date_time_convert = {'TIME': lambda x: str(x).zfill(4),'DATE': lambda x : str(x)}
            self.df = pd.read_csv(fname,delim_whitespace=True,converters=date_time_convert)
            try:
                # Now read in new match,diffs made using lam2/dSI test
                pmatch_l2 = '{}_matches.pairs'.format(fname.split('.')[0])
                pdiff_l2 ='{}_diffs.pairs'.format(fname.split('.')[0])
                pmatch_dsi = '{}_matches_dsi.pairs'.format(fname.split('.')[0])
                pdiff_dsi = '{}_diffs_dsi.pairs'.format(fname.split('.')[0])
                self.matches = pd.read_csv(pmatch_l2,delim_whitespace=True,converters=date_time_convert)
                self.diffs = pd.read_csv(pdiff_l2,delim_whitespace=True,converters=date_time_convert)
                self.matches_dsi = pd.read_csv(pmatch_dsi,delim_whitespace=True,converters=date_time_convert)
                self.diffs_dsi = pd.read_csv(pdiff_dsi,delim_whitespace=True,converters=date_time_convert)
                self.uID = pd.read_csv('{}_uID_l2.pairs'.format(fname.split('.')[0]),delim_whitespace=True,converters=date_time_convert)
                # self.null_split = pd.read_csv('{}_null_split.pairs'.format(fname.split('.')[0]),delim_whitespace=True,converters=date_time_convert)
                self.null_split_diffs = pd.read_csv('{}_diffs_null_split.pairs'.format(fname.split('.')[0]),delim_whitespace=True,converters=date_time_convert)
                self.null_split_matches = pd.read_csv('{}_matches_null_split.pairs'.format(fname.split('.')[0]),delim_whitespace=True,converters=date_time_convert)
                self.nulls = pd.read_csv('{}_nulls.pairs'.format(fname.split('.')[0]),delim_whitespace=True,converters=date_time_convert)
            except FileNotFoundError:
                print('Match or Diff file not found. Are you using synhetics maybe??')
                print('Tried (e.g.,) {}_diffs.pairs'.format(fname.split('.')[0]))

        elif file is False:
            print('Expecting df input')
            self.df = df

    ########## UTILS ###############

    def _surf(self,ax,l2n,b_on_w=True,ylab=True,**kwargs):

        '''
        A function to plot a lambda2 surface over phi,dt (dtmax=4.0). Axis object is returned so other things such as a Title can be added.
        '''
        if b_on_w == True:
            c = 'black'
        else:
            c = 'white'
        params = {
            'savefig.dpi': 500,  # to adjust notebook inline plot size
            'axes.labelsize': 14, # fontsize for x and y labels (was 10)
            'axes.titlesize': 16,
            'axes.edgecolor': c, #white for pres. version
            'font.size': 12, # was 10
            'text.color': c, # change to white for presentations
            'axes.labelcolor' : c, #white for pres. version
            'xtick.color' : c, #white for pres. version
            'ytick.color' : c, #white for pres. version
            'legend.fontsize': 12, # was 10
            'xtick.labelsize': 14,
            'ytick.labelsize': 14,
        }
        matplotlib.rcParams.update(params)
        C  =  ax.contour(self.T,self.F,l2n,levels=[2,5,10,15,20],colors=c)
        Ca =  ax.contour(self.T,self.F,l2n,levels=[1],linewidths=3,label='95% contour',colors=c)
        # Add labels to contours
        ax.clabel(C,C.levels,inline=True,fmt ='%2.0f')
        ax.clabel(Ca,Ca.levels,inline=True,fmt ='%2.0f')
        ax.set_xlim([0, 4.0])
        ax.set_ylim([-90,90])
        ax.set_yticks([-90,-60,-30,0,30,60,90])
        ax.set_xticks([0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0])
        ax.set_xlabel(r'Lag time, $\delta t$ (s)')
        if ylab == True:
            ax.set_ylabel(r'Fast direction, $\phi$, ($\degree$)')

        return

    def _ppm(self,ax,e,n, **kwargs):
        '''
        e - east component seismogram (st[0].data, not stream channel)
        n - north component seismogram (st[1].data, not stream channel)
        '''
        # return C, Ca

        ax.plot(n,e)

        # set limits
        lim = np.abs(max([e.max(),n.max()])) * 1.1
        if 'lims' not in kwargs: kwargs['lims'] = [-lim, lim]
        ax.set_aspect('equal')
        ax.set_ylim(kwargs['lims'])
        ax.set_xlim(kwargs['lims'])

        # set axes labels
        ax.set_xlabel('N')
        ax.set_ylabel('E')

        # turn off tick labels
        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])
        return

    def _mod_tt(self,tr,phase):
        """
        Function to run TauP traveltime models for the SKS phase.
        Returns SKS predictided arrivals (seconds), origin time of the event (t0) as a UTCDateTime obejct and the SKS arrival as a UTCDateTime object
        tr - trace object for which SKS arrival time will be predicted
        """
        model = obspy.taup.tau.TauPyModel(model="iasp91")

        # Add in a test for the case where depth has been gien in meters
        if tr.stats.sac.evdp >= 1000:
            #This might be a bit generous but my events shouldnt be shallower the ~ 1 km or deeper than 1000km anyway (and if this is the case than there is something SERIOUSLY wrong with our earth models)
            traveltime = model.get_travel_times((tr.stats.sac.evdp/1000),tr.stats.sac.gcarc,[phase])[0].time
            print(traveltime)
        elif tr.stats.sac.evdp == 0: # Theres an event where the event data couldnt be found so evdp was set to be 0
            # Having a depth of zero will give us problems so NOW change it to 10.0km exactly (these traveltimes could be very dodgy)
            traveltime = model.get_travel_times(10,tr.stats.sac.gcarc,[phase])[0].time
        else:
            tt = model.get_travel_times((tr.stats.sac.evdp),tr.stats.sac.gcarc,[phase])
            try:
                traveltime = tt[0].time
            except IndexError:
                print('index Error')
                traveltime =None
        return traveltime

    def lam2_surface(self,fstem=None,stat=None,date=None,time=None):
        ''' Function to read  SKS and SKKS .lam2 surface files from sheba
        If syn if False (i.e real data is being used.) Then fstem in needed
        IF syn is True then f1 , f2 are needed
        '''
        # print(fstem)

        t_sks = '{}/SKS/{}_{}_{}??_SKS.lamR'.format(fstem,stat,date,time)
        t_skks = '{}/SKKS/{}_{}_{}??_SKKS.lamR'.format(fstem,stat,date,time)

        sks =glob(t_sks)
        skks = glob(t_skks)
        # print(sks)
        if len(sks) == 0:
            # print('{}/SKS/{}*_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),fstem.split('/')[-1]))
            sks = glob('{}/SKS/{}_{}*_SKS.lamR'.format(fstem,stat,date))
            skks = glob('{}/SKKS/{}_{}*_SKKS.lamR'.format(fstem,stat,date))
            print(sks)

        self.sks_lam2 = np.loadtxt(sks[0])#,skiprows=4) # skip rows not needed for .lamR files
        self.skks_lam2 = np.loadtxt(skks[0])#,skiprows=4)

        self.nfast,self.nlag = self.sks_lam2.shape ;
        self.lag_max = 4.
        [self.T,self.F] = np.meshgrid(np.linspace(0,self.lag_max,num=self.nlag),np.arange(-90,91,1)) ;
        return

    def package(self,outdir,mypath='/Users/ja17375/DiscrePy/Sheba/SAC'):
        '''
        Function to package the SAC files for the observations made into a new directory (for sharing data and results)
        '''
        if os.path.isdir('{}{}'.format(mypath,outdir)) is False:
            #If outdir doesnt exist then Make it ad underlying structure (an additional SAC directory )
            os.makedirs('{}{}/SAC'.format(mypath,outdir))

        def mk_sacfiles(path,stat,date,time):
            '''Function to generate sacfile names'''
            #path='/Users/ja17375/DiscrePy/Sheba/SAC'
            return '{}/{}/{}_{}_{}*BH?.sac'.format(path,stat,stat,date,time)

        def cp_sacfile(sacfile,path,outdir):
            '''Function to copy sacfile to the output directory
            Path - path to the output directorys
            Outdir - output directory name
            '''
            call('cp {} {}/{}/SAC/'.format(sacfile,path,outdir),shell=True)
            print('cp {} {}/{}/SAC/'.format(sacfile,path,outdir))

        mk_sacfiles_mypath = lambda stat,date,time : mk_sacfiles('/Users/ja17375/DiscrePy/Data/SAC_files',stat,date,time)
        cp_sacfiles_mypath = lambda sacfile: cp_sacfile(sacfile,mypath,outdir)
        sacfiles = list(map(mk_sacfiles_mypath,self.P.STAT,self.P.DATE,self.P.TIME))

        for file in sacfiles:
            cp_sacfile(file,mypath,outdir)


    def discrepancy_plot(self,surf_path=None,nplots=2,surfs_to_plot=None,save=False,sigma=1,**kwargs):
        '''Top level plotting function for surfaces to look for discrepancy in Splitting
            surf_path [str] - path to the runs directory containing the surfaces to plot
            nplots - the number of plots that you want (if the surfs_to_plot is not specified)
            surfs_to_plot - allows
            sigma - multiplier to error bounds of splitting results
        '''
        self.spath = surf_path
        self.p_sorted = self.df.sort_values(by='LAM2_BAR',ascending=True)
        self.p_sorted.reset_index(drop=True)
        # Find indicies of events we want to plot
        if surfs_to_plot is None:
            self.surfs= np.round(np.arange(0,len(self.p_sorted),round((len(self.p_sorted)/nplots))))
        else:
            self.surfs = list(set([i if i < len(self.p_sorted) else (len(self.p_sorted)-1) for i in surfs_to_plot])) # This makes sure that indicies are always within the range of available surfaces (stops errors for occuring)
            #self.surfs = self.surfs.sort() # Sorts list in ascending order, has to be done speratly as sort acts of list and returns nothing
        # print(self.surfs)
        if save is True:
            dir = input('Enter Directory you want to save stacked surfaces to > ')

            if os.path.isdir('/Users/ja17375/DiscrePy/Figures/Stacked_Surfaces/{}'.format(dir)) is False:
                make_d = input('Directory {} does not exist, do you want to create it? (y/n)'.format(dir))
                if make_d =='y':
                    print('Ok, creating directory {}'.format(dir))
                    os.mkdir('/Users/ja17375/DiscrePy/Figures/Stacked_Surfaces/{}'.format(dir))
                else:
                    print('Exiting....')
                    sys.exit()

        for s in self.surfs:
            stat,date,time = self.p_sorted.STAT.values[s],self.p_sorted.DATE.values[s], self.p_sorted.TIME.values[s]
            # Note that dtlag and dfast are multiplied through by sigma HERE !
            fast_sks,dfast_sks,lag_sks,dlag_sks = self.p_sorted.FAST_SKS.values[s],(sigma*self.p_sorted.DFAST_SKS.values[s]),self.p_sorted.TLAG_SKS.values[s],(sigma*self.p_sorted.DTLAG_SKS.values[s])
            fast_skks,dfast_skks,lag_skks,dlag_skks = self.p_sorted.FAST_SKKS.values[s],(sigma*self.p_sorted.DFAST_SKKS.values[s]),self.p_sorted.TLAG_SKKS.values[s],(sigma*self.p_sorted.DTLAG_SKKS.values[s])
            lam2,lam2a_sks,lam2a_skks,lam2a_sum = self.p_sorted.LAM2_BAR.values[s],self.p_sorted.LAM2A_SKS.values[s],self.p_sorted.LAM2A_SKKS.values[s],self.p_sorted.LAM2_SUM.values[s]
            l_path = '{}/{}'.format(self.spath,stat) #Path to lambda 2 surfaces for SKS and SKKS
            self.lam2_surface(l_path,stat,date,time)
            fig, (ax0,ax1,ax2) = plt.subplots(1,3,figsize=(24,7),sharey=True)
            fig.patch.set_facecolor('None')
            # ax0.patch.set_facecolor('black')
            if self.syn is True:
                plt.suptitle(r'Syn Stack E1: {} E2: {} $\lambda _2$ value = {:4.3f}'.format(self.syn1[s],self.syn2[s],self.df.LAM2_BAR.values[s]),fontsize=20)
            else:
                plt.suptitle(r'Station {} Date {} Time {}. $\bar{{\lambda_2}}$ = {:4.3f}, $\Delta SI$ = {}'.format(stat,date,time,self.p_sorted.LAM2_BAR.values[s],self.p_sorted.D_SI_Pa.values[s]),fontsize=20)

            self._surf(ax0,self.sks_lam2/lam2a_sks)
            ax0.set_title(r'$\Lambda _2 (\phi,\delta t)$ for SKS')
            self._surf(ax1,self.skks_lam2/lam2a_skks,ylab=False)
            ax1.set_title(r'$\Lambda _2 (\phi,\delta t)$ for SKKS')
            #Plot SKS Solution
            ax0.plot(lag_sks,fast_sks,'.',color='xkcd:blue',label='SKS Solution',markersize=18)
            print('Lag sks {}. Fast SKS {}.'.format(lag_sks,fast_sks))
            ax0.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'-',color='xkcd:blue',linewidth=3)
            ax0.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'-',color='xkcd:blue',linewidth=3)

            #Plot SKKS Solution
            ax0.plot(lag_skks,fast_skks,'.',color='orange',label='SKKS Solution',markersize=18)
            ax0.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'-',color='orange',linewidth=3)
            ax0.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'-',color='orange',linewidth=3)

            #Plot SKS Solution
            ax1.plot(lag_sks,fast_sks,'.',color='xkcd:blue',label='SKS Solution',markersize=18)
            ax1.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'-',color='xkcd:blue',linewidth=3)
            ax1.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'-',color='xkcd:blue',linewidth=3)
            #Plot SKKS Solution
            ax1.plot(lag_skks,fast_skks,'.',color='orange',label='SKKS Solution',markersize=18)
            ax1.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'-',color='orange',linewidth=3)
            ax1.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'-',color='orange',linewidth=3)

            stk = (self.sks_lam2 + self.skks_lam2) / (lam2a_sks + lam2a_skks)
            jf,jt  = np.unravel_index(stk.argmin(),stk.shape)
            stk_fast = np.arange(-90,91,1)[jf]
            stk_lag = np.arange(0,4.025,0.025)[jt]
            self._surf(ax2,stk,ylab=False)


            print('Lam2 BAR is ',stk.min())

            ##########################
            # self.stk_dlag = self.stk_dlag*sigma
            # self.stk_dfast = self.stk_dfast*sigma
            ############################
            ax2.set_ylim([-90,90])
            ax2.set_xlim([0,4])
            ax2.set_yticks([-90,-60,-30,0,30,60,90])
            #Plto SKS solution
            ax2.plot(lag_sks,fast_sks,'.',color='xkcd:blue',label='SKS Solution',markersize=18)
            ax2.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'-',color='xkcd:blue',linewidth=3)
            ax2.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'-',color='xkcd:blue',linewidth=3)
            #Plot SKKS solution
            ax2.plot(lag_skks,fast_skks,'.',color='orange',label='SKKS Solution',markersize=18)
            ax2.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'-',color='orange',linewidth=3)
            ax2.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'-',color='orange',linewidth=3)
            # Plot Stacked solution on SKS
            ax0.plot(stk_lag,stk_fast,'.',label='Stacked Solution',color='xkcd:green',markersize=18)
            # Plot Stacked Solution on SKKS surface
            ax1.plot(stk_lag,stk_fast,'.',label='Stacked Solution',color='xkcd:green',markersize=18)
            ax2.plot(stk_lag,stk_fast,'.',color='xkcd:green',markersize=18)
            ## Add a legend (on ax0)
            # leg = ax0.legend(bbox_to_anchor=(0,1),loc='upper right')
            # plt.setp(leg.get_texts(),color='k')
            ax2.set_title(r'$ \bar{\Lambda_2} (\phi,\delta t)$')
            if save is True:
                # dir = input('Enter Directory you want to save stacked surfaces to > ')
                plt.savefig('/Users/ja17375/DiscrePy/Figures/Stacked_Surfaces/{}/LAM2_{:4.4f}_STAT_{}.eps'.format(dir,lam2,stat),format='eps',dpi=600,transparent=True)
                plt.close()
            elif save is False:
                plt.show()


    def plot_single_phase(self,i,phase,rdir='E_pacific',save=False):
        '''Function to find a SWS result from SHEBA and make a nicer version of the output plot (Lam2 surface, ppm ) '''

        evt = self.df.iloc[i].copy()
        s = '{}/Sheba/Runs/{}/{}/{}/{}_{}_{}??_{}.BH*'.format(self.rpath,rdir,evt.STAT,phase,evt.STAT,evt.DATE,evt.TIME,phase)
        s_c = '{}/Sheba/Runs/{}/{}/{}/{}_{}_{}??_{}_corr.BH*'.format(self.rpath,rdir,evt.STAT,phase,evt.STAT,evt.DATE,evt.TIME,phase)
        st = obspy.read(s) # Read in uncorrected trace (trimmed and output from Sheba.py)
        d_st = st[0].stats.delta # Read the sampling interval
        tt = self._mod_tt(st[0],phase) - 60# Get traveltime (minus 60 seconds as thats when traces start)
        st_c = obspy.read(s_c)
        t = tt + np.arange(0,st[0].stats.npts*d_st,d_st)

        self.lam2_surface('{}/Sheba/Runs/{}/{}'.format(self.rpath,rdir,evt.STAT),evt.STAT,evt.DATE,evt.TIME)

        if phase == 'SKS':
            wbeg = evt.WBEG_SKS
            wend = evt.WEND_SKS
            l2n = self.sks_lam2/evt.LAM2A_SKS
            phi = evt.FAST_SKS
            dphi = evt.DFAST_SKS
            dt = evt.TLAG_SKS
            ddt = evt.DTLAG_SKS
            spol = evt.SPOL_SKS
            Q = evt.Q_SKS

        elif phase == 'SKKS':
            wbeg = evt.WBEG_SKKS
            wend = evt.WEND_SKKS
            l2n = self.skks_lam2/evt.LAM2A_SKKS
            phi = evt.FAST_SKKS
            dphi = evt.DFAST_SKKS
            dt = evt.TLAG_SKKS
            ddt = evt.DTLAG_SKKS
            spol = evt.SPOL_SKKS
            Q = evt.Q_SKKS
        ####### Plotting ########
        fig = plt.figure(figsize=(18,12))
        gs = GridSpec(4,3)
        ax1 = fig.add_subplot(gs[0,:])   # top row for traces
        ax2 = fig.add_subplot(gs[1,:])   # top row for traces
        ax3 = fig.add_subplot(gs[2,0])   # box for uncorrected particle motion
        ax4 = fig.add_subplot(gs[3,0])   # box for corrected particle motion
        ax5 = fig.add_subplot(gs[2:,1:]) # 2x2 box for lam2 surface
        #Plot traces
        ax1.plot(t,st[0].data,'-')
        ax1.plot(t,st[1].data,'-',color='darkorange')
        ax1.axvline(wbeg,linewidth=1.5,color='r')
        ax1.axvline(wend,linewidth=1.5,color='r')
        ax1.set_xlim([wbeg - 20,wend+20])

        # Plot corrected traces
        ax2.plot(t,st_c[0].data,'--')
        ax2.plot(t,st_c[1].data,'--')
        ax2.axvline(wbeg,linewidth=1.5,color='r')
        ax2.axvline(wend,linewidth=1.5,color='r')
        ax2.set_xlim([wbeg - 20,wend + 20])

        #Plot Uncorrected particle motion
        st2 = st.copy()
        st2[0].trim(st2[0].stats.starttime + (wbeg - tt), st2[0].stats.starttime + (wend - tt))
        st2[1].trim(st2[1].stats.starttime + (wbeg - tt), st2[1].stats.starttime + (wend - tt))
        self._ppm(ax3,st2[0].data,st2[1].data)
        # Plor corrected particle motion
        st2_c = st.copy()
        st2_c[0].trim(st2_c[0].stats.starttime + (wbeg - tt), st2_c[0].stats.starttime + (wend - tt))
        st2_c[1].trim(st2_c[1].stats.starttime + (wbeg - tt), st2_c[1].stats.starttime + (wend - tt))
        self._ppm(ax4,st2_c[0].data,st2_c[1].data)
        #

        self._surf(ax5,l2n)
        ax5.plot(dt,phi,'X',markersize=10,color='r')

        #### Just plot the Lam2 surface
        fig2 = plt.figure(figsize=(9,9))
        ax = fig2.add_subplot(111)

        y = st[0].stats.starttime.year
        d =st[0].stats.starttime.julday
        h = st[0].stats.starttime.hour
        m = st[0].stats.starttime.minute
        s = st[0].stats.starttime.second

        self._surf(ax,l2n)
        ax.plot(dt,phi,'X',markersize=16,color='r')
        plt.suptitle(r'Stat: {}, Date: {:04}{:03} Time: {:02}:{:02}:{:02}'.format(evt.STAT,y,d,h,m,s),fontsize=14)
        plt.title(r'$\phi$ = {:2.2f}$\degree$ +/- {:2.2f}, $\delta t$ = {:2.2f}s +/- {:2.2f}, SPOL = {:3.1f}$\degree$, Q = {:2.3f}'.format(phi,dphi,dt,ddt,spol,Q),fontsize=14)

        if save == True:
            fig2.savefig('{}Figures/{}_null_ex.eps'.format(self.rpath,evt.STAT),dpi = 600)

        plt.show()
        # return st

    def plot_dist_v_discrep(self):

        fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize = [12,6])

        ax1.plot(self.matches.DIST,self.matches.LAM2_BAR,marker='.',color='blue',label='Matching')
        ax1.plot(self.diffs.DIST,self.diffs.LAM2_BAR,marker='.',color='darkorange',label='Discrepant')
        ax1.set_xlim([105, 140])
        ax1.set_ylim([0, np.around(self.df.LAM2_BAR.max(),decimals=1)])
        ax1.set_xlabel('Epicentral Distance (Deg)')
        ax1.set_ylabel(r'$\bar{\lambda _2}$ values')
        ax1.set_title('Match/Diff according to  2 sigma test')
        ax1.legend()

        ax2.plot(self.matches.DIST,self.matches.D_SI_Pr,marker='.',color='blue',label='Matching')
        ax2.plot(self.diffs.DIST,self.diffs.D_SI_Pr,marker='.',color='darkorange',label='Discrepant')
        ax2.set_xlim([105,140])
        ax2.set_ylim([0,np.around(self.df.D_SI_Pr.max(),decimals=1)])
        ax2.set_xlabel('Epicentral Distrance (Deg)')
        ax2.set_ylabel(R'$\Delta SI $')

        ax3.plot(self.matches.DIST,self.matches.LAM2_BAR,marker='.',color='blue',label='Matching')
        ax3.plot(self.diffs.DIST,self.diffs.LAM2_BAR,marker='.',color='darkorange',label='Discrepant')
        ax3.set_xlim([105, 140])
        ax3.set_ylim([0, np.around(self.df.LAM2_BAR.max(),decimals=1)])
        ax3.set_xlabel('Epicentral Distance (Deg)')
        ax3.set_ylabel(r'$\bar{\lambda _2}$ values')
        ax3.set_title(r'Match/Diff according to $\bar{lambda _2} / \Delta SI $ test')
        ax3.legend()

        ax4.plot(self.matches.DIST,self.matches.D_SI_Pr,marker='.',color='blue',label='Matching')
        ax4.plot(self.diffs.DIST,self.diffs.D_SI_Pr,marker='.',color='darkorange',label='Discrepant')
        ax4.set_xlim([105,140])
        ax4.set_ylim([0,np.around(self.df.D_SI_Pr.max(),decimals=1)])
        ax4.set_xlabel('Epicentral Distrance (Deg)')
        ax4.set_ylabel(R'$\Delta SI $')

        plt.show()

    def spol_v_baz(self,save=False,fname=None):
        '''Plot measure spource polarisation against back azimuth (a proxy for source pol). In theory SPOL ~= BAZ  '''

        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(16,8))

        ax1.plot(self.df.SPOL_SKS,self.df.BAZ,marker='.',color='blue')
        ax1.plot(self.df.SPOL_SKKS,self.df.BAZ,marker='x',color='blue')
        # ax1.set_xlim()
        ax2.plot(self.diffs.SPOL_SKS,self.diffs.BAZ,marker='.',color='darkorange',linestyle='None')
        ax2.plot(self.matches.SPOL_SKKS,self.matches.BAZ,marker='x',color='darkorange')


        plt.show()

    def l2_v_SNR(self,fname=None,save=False):
        '''Plots lambda 2 for clear split pairs (Q> 0.5 for both phases) against the pair SNR. '''
        fig =  plt.figure(figsize=(18,12))
        # Isolate clear split pairs

        splits = self.df[(self.df.Q_SKS >= 0.5) & self.df.Q_SKKS >= 0.5]
        # Calc snr
        SNR = (splits.SNR_SKS + splits.SNR_SKKS )/2
        # Define a range of SNR values to evaltue models for ...
        snr_mod = np.arange(SNR.min()-0.1,SNR.max()+20,0.1)
        #print(len(snr_mod))
        # plot lam2 v SNR
        ax1 = fig.add_subplot(121)
        ax1.scatter(SNR,splits.LAM2_BAR,marker='.')
        ax1.set_title(r'$\bar{\lambda_2}$ for split pairs against SNR')
        ax1.set_xlabel('Signal-to-Noise Ratio (SNR)')
        ax1.set_ylabel(r'$\bar{\lambda_2}$')
        #plot log(lam2) v SNR
        # Semiplog plot removed as it wasnt determined to be unhelpful
        ax2 = fig.add_subplot(122)
        ax2.loglog(SNR,splits.LAM2_BAR,'.')
        # Fit a straight line to THIS (log-log) data
        m2, c2 = np.polyfit(np.log(SNR),np.log(splits.LAM2_BAR),1)
        # Calculate line as estimated by ployfit
        y_fit2 = m2 * np.log(snr_mod) + (c2)
        # Also calulate the model we got from the syntetics (m = ~2, c = 0.37)
        m = -2.04
        c = 0.37
        y_mod_log = (m * np.log(snr_mod)) + c
        print(max(SNR),min(SNR))

        ax2.set_title(r'log ($\bar{\lambda_2})$ for split pairs against log(SNR)')
        ax2.set_xlabel('Signal-to-Noise Ratio (SNR)')
        ax2.set_ylabel(r'$(\bar{\lambda_2})$')
        ax2.text(0.06,0.2,r'Data Model : $ln(\lambda _2$) = {:4.2f} * ln(SNR) + {:4.2f}'.format(m2,c2),transform=ax2.transAxes)
        ax2.text(0.06,0.15,r'Synthetics Model : $ln(\lambda _2$) = -2.04 * ln(SNR) + 0.37'.format(m2,c2),transform=ax2.transAxes)

        #Add Fit lines to linear plot

        f1, = ax1.plot(snr_mod,np.exp(y_fit2),'k--',label='Log fit (to data)')
        f2, = ax1.plot(snr_mod,np.exp(y_mod_log),'r--',label='Model from Synthetics')

        ax2.loglog(snr_mod,np.exp(y_mod_log),'--r')
        ax2.loglog(snr_mod,np.exp(y_fit2),'--k')
        ax1.legend((f1,f2),('Model from Data','Model from Synthetics'))

        if save is True:
            plt.savefig('/Users/ja17375/DiscrePy/Figures/{}.eps'.format(fname),format='eps',dpi=1000)
        plt.show()

    def l2_dSI_SNR(self,fname=None,save=False):
        '''Plots lam2 and dSI agaisnt each other for SKS and SKKS, coloured by signal-to-noise ratio (SNR)'''
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        # Isolate clear split pairs
        splits = self.df[(self.df.Q_SKS >= 0.5) & self.df.Q_SKKS >= 0.5]
        print(len(splits))
        c1 = ax1.scatter(splits.LAM2_BAR,splits.D_SI_Pr,c=splits.SNR_SKS,marker='.',vmax=15,vmin=5)
        plt.colorbar(c1,ax=ax1)
        c2 = ax2.scatter(splits.LAM2_BAR,splits.D_SI_Pr,c=splits.SNR_SKKS,marker='.',vmax=15,vmin=5)
        plt.colorbar(c2,ax=ax2)
        #Plot Thresholds
        ax1.plot([0,max([np.around(splits.LAM2_BAR.max(),decimals=2),np.around(splits.LAM2_BAR.max(),decimals=2) ] ) ],[0.2,0.2],'k--')
        ax1.plot([0.03,0.03],[0,max([np.around(splits.D_SI_Pr.max(),decimals=2),np.around(splits.D_SI_Pr.max(),decimals=2)])],'k--')
        ax1.set_ylim([0,max([np.around(splits.D_SI_Pr.max(),decimals=1),np.around(splits.D_SI_Pr.max(),decimals=2)])])
        ax1.set_xlim([0,max([np.around(splits.LAM2_BAR.max(),decimals=1),np.around(splits.LAM2_BAR.max(),decimals=2)])])
        ax1.set_xlabel(r'$\bar{\lambda _2}$')
        ax1.set_ylabel(r'$\Delta$ SI')
        ax1.set_title('SKS')
        ax2.plot([0,max([np.around(splits.LAM2_BAR.max(),decimals=2),np.around(splits.LAM2_BAR.max(),decimals=2) ] ) ],[0.2,0.2],'k--')
        ax2.plot([0.03,0.03],[0,max([np.around(splits.D_SI_Pr.max(),decimals=2),np.around(splits.D_SI_Pr.max(),decimals=2)])],'k--')
        ax2.set_ylim([0,np.around(splits.D_SI_Pr.max(),decimals=2)])
        ax2.set_xlim([0,np.around(splits.LAM2_BAR.max(),decimals=2)])
        ax2.set_xlabel(r'$\bar{\lambda _2}$')
        ax2.set_ylabel(r'$\Delta$ SI')
        ax2.set_title('SKKS')
        if save is True:
            plt.savefig('/Users/ja17375/DiscrePy/Figures/{}.eps'.format(fname),format='eps',dpi=1000)
        plt.show()

    def l2_dSI_hist(self,save=False):

        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,6))
        bins_l2 = np.linspace(0,0.5,20)
        bins_dsi = np.linspace(0,3.5,20)
        m_splits = self.matches[(self.matches.Q_SKS > 0.5) & (self.matches.Q_SKKS > 0.5)]
        d_splits = self.diffs[(self.diffs.Q_SKS > 0.5) & (self.diffs.Q_SKKS > 0.5)]

        ax1.hist([m_splits.LAM2_BAR,d_splits.LAM2_BAR],bins_l2, histtype='bar', stacked=True,label=["'Matching'","'Discrepent'"])
        ax1.set_xlabel(r'$\bar {\lambda _2}$ values')
        ax1.set_ylabel('Count')
        ax1.legend()
        ax1.set_xlim([0, np.around(self.df.LAM2_BAR.max(),decimals=1)])
        ax2.hist([m_splits.D_SI_Pr,d_splits.D_SI_Pr],bins_dsi, histtype='bar', stacked=True,label=["'Matching'","'Discrepent'"])
        ax2.set_xlabel(r'$\Delta SI$ values')
        ax2.set_ylabel('Count')
        ax2.legend()
        ax2.set_xlim([0,np.around(self.df.D_SI_Pr.max(),decimals=1)])
        # ax2.set_ylim([0,4.0])
        if save is True:
            plt.savefig('/Users/ja17375/DiscrePy/Figures/hist_and_dVs.png',format='png',dpi=1000)

        plt.show()

    def plot_dist_v_split(self,save=False):

        # Plot colored by Lambda 2
        fig1,(ax1,ax2) = plt.subplots(1,2, figsize = [12,6])

        deltaF = 90 -abs(abs(self.df.FAST_SKS - self.df.FAST_SKKS) - 90 ) # This should allow for differences to vary between 0 and 90 whilst also dealing with the wrap around properly
        ddt = abs(self.df.TLAG_SKS) - abs(self.df.TLAG_SKKS)

        C1 = ax1.scatter(self.df.DIST,deltaF,c=self.df.LAM2_BAR,marker='.',label=r'$| \phi_{SKS} - \phi_{SKKS} |')

        ax1.set_xlim([105, 140])
        ax1.set_ylim([0,90])
        ax1.set_xlabel('Epicentral Distance (Deg)')
        ax1.set_ylabel(r'$ |\phi_{SKS} - \phi_{SKKS}|$')
        ax1.set_title(r'Difference in $\phi$, coloured by $\Lambda$')

        C2 = ax2.scatter(self.df.DIST,ddt,c=self.df.LAM2_BAR,marker='.',label='delta deltaT')
        #ax2.plot(df.diffs.DIST,df.diffs.D_SI,'r.',label='Discrepant')
        ax2.set_xlim([105,140])
        ax2.set_ylim([0,4])
        ax2.set_xlabel('Epicentral Distrance (Deg)')
        ax2.set_ylabel(r'$| \delta t_{SKS} - \delta t_{SKKS} | $')
        ax2.set_title(r'Difference in $(\delta_t)$, coloured by $\Lambda$')
        cbar1 = fig1.colorbar(C1,use_gridspec=True)
        cbar1.set_label(r'$\Lambda$',rotation=270)

        fig2,(ax3) = plt.subplots(figsize= [6,6])
        C3 = ax3.scatter(ddt,(deltaF),c=self.df.LAM2_BAR,marker='.')
        ax3.set_xlim([0,4])
        ax3.set_ylim([0, 90])
        ax3.set_ylabel(r'$| \phi_{SKS} - \phi_{SKKS} | $')
        ax3.set_xlabel(r'$| \delta t_{SKS} - \delta t_{SKKS} | $')
        ax3.set_title(r'Difference in $\phi$ v difference in $\delta t$ coloured by $\Lambda$')

        cbar2 = fig2.colorbar(C3,use_gridspec=True)
        cbar2.set_label(r'$\Lambda$',rotation=270)
        # plot coloured by D_SI
        fig3,(ax4,ax5) = plt.subplots(1,2, figsize = [14,6])

        C4 = ax4.scatter(self.df.DIST,deltaF,c=self.df.D_SI,cmap='magma',marker='.',label=r'$| \phi_{SKS} - \phi_{SKKS} |')

        ax4.set_xlim([105, 140])
        ax4.set_ylim([0,90])
        ax4.set_xlabel('Epicentral Distance (Deg)')
        ax4.set_ylabel(r'$ |\phi_{SKS} - \phi_{SKKS}|$')
        ax4.set_title(r'Difference in $\phi$, coloured by $\Delta SI$')
        C5 = ax5.scatter(self.df.DIST,ddt,c=self.df.D_SI,cmap='magma',marker='.',label='delta deltaT')
        #ax2.plot(df.diffs.DIST,df.diffs.D_SI,'r.',label='Discrepant')
        ax5.set_xlim([105,140])
        ax5.set_ylim([0,4])
        ax5.set_xlabel('Epicentral Distrance (Deg)')
        ax5.set_ylabel(r'$| \delta t_{SKS} - \delta t_{SKKS} | $')
        ax5.set_title(r'Difference in $\delta_t$, coloured by $\Delta SI$')
        cbar3 = fig3.colorbar(C5,use_gridspec=True)
        cbar3.set_label(r'$\Delta SI$',rotation=270)

        fig4,(ax6) = plt.subplots(figsize= [6,6])
        C6 = ax6.scatter(ddt,(deltaF),cmap='magma',c=self.df.D_SI,marker='.')
        ax6.set_xlim([0,4])
        ax6.set_ylim([0, 90])
        ax6.set_ylabel(r'$| \phi_{SKS} - \phi_{SKKS} | $')
        ax6.set_xlabel(r'$| \delta t_{SKS} - \delta t_{SKKS} | $')
        ax6.set_title(r'Difference in $\phi$ v difference in $\delta t$ coloured by $\Delta SI$')
        cbar4 = fig4.colorbar(C6,use_gridspec=True)
        cbar4.set_label(r'$\Delta SI$',rotation=270)

        if save is True:
            fig1.savefig('/Users/ja17375/DiscrePy/Sheba/Figures/diff_v_dist_Lambda.eps',format='eps',transparent=True,dpi=400)
            fig2.savefig('/Users/ja17375/DiscrePy/Sheba/Figures/delfast_v_deltlag_Lambda.eps',format='eps',transparent=True,dpi=400)
            fig3.savefig('/Users/ja17375/DiscrePy/Sheba/Figures/diff_v_dist_DSI.eps',format='eps',transparent=True,dpi=400)
            fig4.savefig('/Users/ja17375/DiscrePy/Sheba/Figures/delfast_v_deltlag_DSI.eps',format='eps',transparent=True,dpi=400)

            plt.close('all')
        else:
            plt.show()

    def Q_v_l2_dSI(self):
        '''
        Make plots of Q factor (Wuestefeld et al) sgainst lambda 2 and Delta SI
        '''
        fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(14,6),sharex=True)#,sharey=True)

        # C1 = ax1.scatter(self.df.Q_SKS,self.df.LAM2,c=self.df.SNR_SKS,marker='.',vmin=2,vmax=20)
        ax1.plot(self.matches.Q_SKS,self.matches.LAM2_BAR,marker='.',color='blue',ls='None')
        ax1.plot(self.diffs.Q_SKS,self.diffs.LAM2_BAR,marker='.',color='darkorange',ls='None')
        ax1.set_xlabel('SKS Q factor')
        ax1.set_ylabel(r'$\bar{\lambda _2}$')
        ax1.set_xlim([-1,1])
        ax1.set_ylim([0,np.around(self.df.LAM2_BAR.max(),decimals=1)])
        # plt.colorbar(C1,ax=ax1)

        ax2.plot(self.matches.Q_SKKS,self.matches.LAM2_BAR,marker='.',color='blue',ls='None')
        ax2.plot(self.diffs.Q_SKKS,self.diffs.LAM2_BAR,marker='.',color='darkorange',ls='None')
        #ax2.plot(self.diffs.Q_SKKS,self.diffs.LAM2,marker='.',color='darkorange',ls='None')
        ax2.set_xlabel('SKKS Q factor')
        ax2.set_xlim([-1,1])
        ax2.set_ylim([0,np.around(self.df.LAM2_BAR.max(),decimals=1)])

        ax3.plot(self.matches.Q_SKS,self.matches.D_SI,marker='.',color='blue',ls='None')
        ax3.plot(self.diffs.Q_SKS,self.diffs.D_SI,marker='.',color='darkorange',ls='None')
        #ax3.plot(self.diffs.Q_SKS,self.diffs.D_SI,marker='.',color='darkorange',ls='None')
        ax3.set_xlabel('SKS Q factor')
        ax3.set_ylabel(r'$\Delta SI$')
        ax3.set_xlim([-1,1])
        ax3.set_ylim([0,np.around(self.df.D_SI.max(),decimals=1)])

        ax4.plot(self.matches.Q_SKKS,self.matches.D_SI,marker='.',color='blue',ls='None')
        ax4.plot(self.diffs.Q_SKKS,self.diffs.D_SI,marker='.',color='darkorange',ls='None')
        #ax4.plot(self.diffs.Q_SKKS,self.diffs.D_SI,marker='.',color='darkorange',ls='None')
        ax4.set_xlabel('SKKS Q factor')
        ax4.set_xlim([-1,1])
        ax4.set_ylim([0,np.around(self.df.D_SI.max(),decimals=1)])
        plt.show()

    def plot_SNR(self):
        '''
        Make plots of SNR v dfast. In the style of phi_i v SNR from Restivo and Helffrich (2006)
        '''
        fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,sharey='row',figsize=(6,6))
        #Plot SNR for SKS
        ax1.plot(self.df.SNR_SKS,self.df.DFAST_SKS,'k.')
        ax1.set_ylabel('dfast')
        ax1.set_xlabel('S/N ratio')
        ax1.set_title('d fast SKS determination dependance on S/N')
        #Plot SNR for SKKS
        ax2.plot(self.df.SNR_SKKS,self.df.DFAST_SKKS,'k.')
        ax2.set_ylabel('dfast')
        ax2.set_xlabel('S/N ratio')
        ax2.set_title('d fast SKKS determination dependance on S/N')

        plt.show()

    def hist_SNR(self):
        ''' Plots SNR histogram'''
        print('Max S/N SKS: ', self.df.SNR_SKS.max())
        print('Max S/N SKKS: ', self.df.SNR_SKKS.max())
        fig, (ax1,ax2) = plt.subplots(1,2,figsize=(14,7))
        ax1.hist(self.df.SNR_SKS,bins=np.arange(0,100,5),histtype='bar')
        ax1.set_xlabel('S/N ratio (SKS)')
        ax1.set_ylabel('Count')
        ax1.set_xlim([0, 100])
        ax2.hist(self.df.SNR_SKKS,bins=np.arange(0,100,5),histtype='bar')
        ax2.set_xlabel('S/N ratio (SKKS)')
        ax2.set_ylabel('Count')
        ax2.set_xlim([0, 100])
        plt.show()

    def plot_SNR_v_l2(self):
        '''
        Make plots of SNR v lambda2. In the style of phi_i v SNR from Restivo and Helffrich (2006)
        '''
        fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,sharey='row',figsize=(8,8))
        #Plot SNR for SKS
        ax1.plot(self.df.SNR_SKS,self.df.LAM2_BAR,'k.')
        ax1.set_ylabel(r'$\bar{\lambda_2}$')
        ax1.set_xlabel('S/N ratio')
        ax1.set_title(r'$\bar{\lambda_2}$ v SNR for SKS')
        #Plot SNR for SKKS
        ax2.plot(self.df.SNR_SKKS,self.df.LAM2_BAR,'k.')
        ax2.set_ylabel(r'$\bar{\lambda_2}$')
        ax2.set_xlabel('S/N ratio')
        ax2.set_title(r'$\bar{\lambda_2}$ v SNR for SKKS')
        plt.show()

    def plot_l2sum_v_l2bar(self,save=False):
        '''
        Function to plot lambda2 bar against the sum of lambda2 for each phase in the pair

        We only plto null split and split pairs
        '''
        fig,ax = plt.subplots(1,1,figsize=(7,7))
        splits = self.diffs.append(self.matches)
        m_splits = self.matches[(self.matches.Q_SKS > 0.5) & (self.matches.Q_SKKS > 0.5)]
        d_splits = self.diffs[(self.diffs.Q_SKS > 0.5) & (self.diffs.Q_SKKS > 0.5)]
        #Plot split pairs
        # ax.scatter((m_splits.LAM2_SKS + m_splits.LAM2_SKKS),m_splits.LAM2_BAR,marker='.',label=r'$2 \sigma$ matching')
        # ax.scatter((d_splits.LAM2_SKS + d_splits.LAM2_SKKS),d_splits.LAM2_BAR,marker='.',label=r'$2 \sigma$ discrepant')
        # ax.scatter((self.df.LAM2_SKS + self.df.LAM2_SKKS),self.df.LAM2_BAR,marker='.',label='All pairs')
        # ax.scatter((splits.LAM2_SKS + splits.LAM2_SKKS),splits.LAM2_BAR,marker='.',label='Split Pairs')
        ax.scatter((self.matches.LAM2_SUM),self.matches.LAM2_BAR,marker='.',label=r'$\bar{\lambda_2} $ Matching')
        ax.scatter((self.diffs.LAM2_SUM),self.diffs.LAM2_BAR,marker='x',label=r'$\bar{\lambda_2} $ Discrepant')
        ax.scatter((self.null_split_matches.LAM2_SUM),self.null_split_matches.LAM2_BAR,marker='.',label=r'$\bar{\lambda_2} $ Matching Null-Split')
        ax.scatter((self.null_split_diffs.LAM2_SUM),self.null_split_diffs.LAM2_BAR,marker='x',label=r'$\bar{\lambda_2} $ Discrepant Null-Split')
        #Plot null - split pairs
        # ax.scatter((self.null_split.LAM2_SKS + self.null_split.LAM2_SKKS),self.null_split.LAM2_BAR,marker='.',c='red',label='Null-Split Pairs')
        mod = np.linspace(0,0.2,10)
        # ax.plot(mod,mod,'k--',label=r'$\lambda_2^{SKS} + \lambda_2^{SKKS} = \bar{\lambda_2}$')
        ax.plot(mod,mod,'k-.',label=r'$\bar{\lambda_2} = (\lambda_2^{SKS 95\%} + \lambda_2^{SKKS 95\%}) $')
        ax.set_xlabel(r'$\lambda_2^{SKS} + \lambda_2^{SKKS}$')
        ax.set_ylabel(r'$\bar{\lambda_2}$')
        ax.set_xlim([0,0.07])
        ax.set_ylim([0,0.1])
        ax.legend(fontsize='medium')
        #plt.colorbar(C)\
        if save == True:
            plt.savefig('/Users/ja17375/Thesis/Lambda2_Paper/Figs/Lam2bar_v_Lam2sum_l2_dsi.png',format='eps',transparent=True,dpi=400) # /Users/ja17375/DiscrePy/Figures/
            plt.show()
        else:
            plt.show()

    def plot_l2_v_dsi(self,save=False):
        '''
        Function to plot lambda2 bar against dsi for pair in the E pacific

        We only plto null split and split pairs
        '''
        # Find pairs in the region of interest (Longitdinally anyway)
        d = self.diffs[(self.diffs.SKS_PP_LON > -170) & (self.diffs.SKS_PP_LON < -80) & (self.diffs.SKKS_PP_LON > -170) & (self.diffs.SKS_PP_LON < -80)]
        m = self.matches[(self.matches.SKS_PP_LON > -170) & (self.matches.SKS_PP_LON < -80) & (self.matches.SKKS_PP_LON > -170) & (self.matches.SKS_PP_LON < -80)]

        fig,ax = plt.subplots(1,1,figsize=(7,7))
        ax.scatter(m.LAM2_BAR,m.D_SI_Pr,marker='.',label=r'$\bar{\lambda_2} $ Matching')
        ax.scatter(d.LAM2_BAR,d.D_SI_Pr,marker='x',label=r'$\bar{\lambda_2} $ Discrepant')
        #Plot null - split pairs
        # ax.scatter((self.null_split.LAM2_SKS + self.null_split.LAM2_SKKS),self.null_split.LAM2_BAR,marker='.',c='red',label='Null-Split Pairs')

        ax.plot([0, 0.2],[0.4,0.4],'k-.',label=r'$\Delta SI$ threshold from Deng et al., (2017)')
        ax.set_xlabel(r'$\bar{\lambda_2}$')
        ax.set_ylabel(r'$\Delta SI$')
        ax.set_title('For the Eastern Pacific Region')
        ax.set_xlim([0.01,np.around(max(d.LAM2_BAR.max(),m.LAM2_BAR.max()), decimals=3)])
        # ax.set_ylim([0,0.2])
        ax.legend(fontsize='medium')
        #plt.colorbar(C)\
        if save == True:
            # plt.savefig('/Users/ja17375/Thesis/Lambda2_Paper/Figs/Lam2bar_v_Lam2sum_l2_only.png',format='png',transparent=True,dpi=400) # /Users/ja17375/DiscrePy/Figures/
            plt.show()
        else:
            plt.show()

    def l2_v_dSI(self):
        '''Plot lambda2 (bar) agaisnt splitting intensity'''
        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(12,6))
        ax1.plot(self.diffs_dsi.LAM2_BAR,self.diffs_dsi.D_SI_Pr,'k.',label='Discrepant')
        ax1.plot(self.matches_dsi.LAM2_BAR,self.matches_dsi.D_SI_Pr,'kx',label='Matching')
        ax2.plot(self.diffs.LAM2_BAR,self.diffs.D_SI_Pr,'k.',label='Discrepant')
        ax2.plot(self.matches.LAM2_BAR,self.matches.D_SI_Pr,'kx',label='Matching')
        ax1.set_xlabel(r'$\bar{\lambda_2} $ value', fontsize=14)
        ax1.set_ylabel(r'$\Delta SI$ value',fontsize=14)
        ax2.set_xlabel(r'$\bar{\lambda_2} $ value', fontsize=14)
        ax2.set_ylabel(r'$\Delta SI$ value',fontsize=14)
        ax1.set_ylim([0,1])
        ax1.set_xlim([0,0.2])
        ax2.set_ylim([0,1])
        ax2.set_xlim([0,0.2])
        ax1.legend(fontsize='medium',loc='best')
        ax2.legend(fontsize='medium',loc='best')
        ax1.set_title('Classified using Deng et al., (2017). [$\Delta SI > 0.4$]')
        ax2.set_title(r'Classified using $\bar{\lambda_2}$ & $\Delta SI > 0.3 $')
        plt.show()

    def phi_dt_diff_latlon(self):
        '''
        Plot the difference between phi and delta t for SKS-SKKS pairs, plotted against the pair midpoint postion (defined by sqrt(LAT^2 + LONG^2))
        '''
        fig,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(8,12))
        # Set up Regional Section of pairs
        d = self.diffs[(self.diffs.SKS_PP_LON > -160) & (self.diffs.SKS_PP_LON < -80)]
        ns = self.null_split[(self.null_split.SKS_PP_LON > -160) & (self.null_split.SKS_PP_LON < -80)]
        m = self.matches[(self.matches.SKS_PP_LON > -160) & (self.matches.SKS_PP_LON < -80)]

        latlon_diff = np.sqrt(((d.SKS_PP_LON.values+d.SKKS_PP_LON.values)/2)**2 + ((d.SKS_PP_LAT.values+d.SKKS_PP_LAT.values)/2)**2)
        latlon_match = np.sqrt(((m.SKS_PP_LON.values+m.SKKS_PP_LON.values)/2)**2 + ((m.SKS_PP_LAT.values+m.SKKS_PP_LAT.values)/2)**2)

        ax1.plot(latlon_diff,abs(abs(d.FAST_SKS)-abs(d.FAST_SKKS)),'k.',label='Discrepant')
        ax2.plot(latlon_diff,abs(d.TLAG_SKS-d.TLAG_SKKS) ,'k.')

        ax1.plot(latlon_match,abs(abs(m.FAST_SKS)-abs(m.FAST_SKKS)),'b.',label='Matching')
        ax2.plot(latlon_match,abs(m.TLAG_SKS-m.TLAG_SKKS) ,'b.')
        #Add legend
        ax1.legend()
        # Add x,y labels
        ax2.set_xlabel(r'(Lat + Long)$^{\frac{1}{2}}$ ')
        ax1.set_ylabel(r'$|\phi_{SKS} - \phi_{SKKS}|$ ($\degree$)')
        ax2.set_ylabel('$| \delta t_{SKS} - \delta t_{SKKS}|$ (s)')
        #set limits
        ax1.set_xlim([80, 170])
        ax1.set_ylim([0, 35])
        ax2.set_ylim([0, 1.6])
        ax1.set_title('Fast Direction and Lag time for Pairs in E. Pac.')
        plt.savefig('/Users/ja17375/DiscrePy/Figures/Phi_dt_v_Long_Lat_sq.png',dpi=400)
        plt.show()

    def phi_dt_diff_lat(self):
        '''
        Plot the difference between phi and delta t for SKS-SKKS pairs, plotted against the pair midpoint postion (defined by sqrt(LAT^2 + LONG^2))
        '''
        fig,((ax1,ax3),(ax2,ax4)) = plt.subplots(2,2,sharey=True,figsize=(8,12))
        # Set up Regional Section of pairs
        d = self.diffs[(self.diffs.SKS_PP_LON < -130)]
        ns_sks = self.null_split[(self.null_split.SKS_PP_LON < -130)  & (self.null_split.Q_SKS > 0.5)]
        ns_skks = self.null_split[(self.null_split.SKS_PP_LON < -130)  & (self.null_split.Q_SKKS > 0.5)]
        m = self.matches[(self.matches.SKS_PP_LON < -130) ]

        lat_diff = (d.SKS_PP_LAT + d.SKKS_PP_LAT)/2
        lat_match = (m.SKS_PP_LAT + m.SKKS_PP_LAT)/2
        lat_ns = (ns_skks.SKS_PP_LAT + ns_skks.SKKS_PP_LAT)/2

        ax1.plot(abs(abs(d.FAST_SKS)-abs(d.FAST_SKKS)),lat_diff,'k.',label='Discrepant')
        ax2.plot(abs(d.TLAG_SKS-d.TLAG_SKKS),lat_diff,'k.')
        ax3.plot(d.FAST_SKS,lat_diff,'g.',label='SKS')
        ax3.plot(d.FAST_SKKS,lat_diff,'gd',label='SKKS')
        ax4.plot(d.TLAG_SKS,lat_diff,'g.',label='SKS')
        ax4.plot(d.TLAG_SKKS,lat_diff,'gd',label='SKKS')

        ax1.plot(abs(abs(m.FAST_SKS)-abs(m.FAST_SKKS)),lat_match,'b.',label='Matching')
        ax2.plot(abs(m.TLAG_SKS-m.TLAG_SKKS),lat_match,'b.')
        ax3.plot(m.FAST_SKS,lat_match,'b.',label='SKS')
        ax3.plot(m.FAST_SKKS,lat_match,'bd',label='SKKS')
        ax4.plot(m.TLAG_SKS,lat_match,'b.',label='SKS')
        ax4.plot(m.TLAG_SKKS,lat_match,'bd',label='SKKS')

        ax1.plot(abs(abs(ns_skks.FAST_SKS)-abs(ns_skks.FAST_SKKS)),lat_ns,'r.',label='Null-Split')
        ax2.plot(abs(ns_skks.TLAG_SKS-ns_skks.TLAG_SKKS),lat_ns,'r.')
        # ax3.plot(lat_ns,ns_skks.FAST_SKS,'r.',label='SKS')
        ax3.plot(ns_skks.FAST_SKKS,lat_ns,'rd',label='SKKS')
        # ax4.plot(lat_ns,ns_skks.TLAG_SKS,'r.',label='SKS')
        ax4.plot(ns_skks.TLAG_SKKS,lat_ns,'rd',label='SKKS')
        #Add legend
        ax1.legend()
        # Add x,y labels
        ax1.set_ylabel(r'Latitude ($\degree$)')
        ax2.set_ylabel(r'Latitude ($\degree$)')
        ax1.set_xlabel(r'$|\phi_{SKS} - \phi_{SKKS}|$ ($\degree$)')
        ax2.set_xlabel('$| \delta t_{SKS} - \delta t_{SKKS}|$ (s)')
        #set limits
        # ax1.set_xlim([80, 170])
        # ax1.set_ylim([0, 35])
        # ax2.set_ylim([0, 1.6])
        ax3.set_xlim([-90,90])
        ax4.set_xlim([0,4.0])
        ax1.set_title('Fast Direction and Lag time for Pairs in E. Pac.')
        plt.savefig('/Users/ja17375/DiscrePy/Figures/Phi_dt_v_Lat.png',dpi=400)
        plt.show()

    def DSI_LAM2_lat(self):
        '''
        Plot the difference between phi and delta t for SKS-SKKS pairs, plotted against the pair midpoint postion (defined by sqrt(LAT^2 + LONG^2))
        '''
        fig,(ax1,ax2) = plt.subplots(2,1,sharey=True,figsize=(8,12))
        params = {
            'savefig.dpi': 400,  # to adjust notebook inline plot size
            'axes.labelsize': 26, # fontsize for x and y labels (was 10)
            'axes.titlesize': 14,
            'font.size': 20, # was 10
            'legend.fontsize': 18, # was 10
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           }

        matplotlib.rcParams.update(params)

        # Set up Regional Section of pairs
        d1 = self.diffs[(self.diffs.SKS_PP_LON < -125) & (self.diffs.SKS_PP_LON >= -140) ]
        ns1 = self.null_split[(self.null_split.SKS_PP_LON < -125) &  (self.null_split.SKS_PP_LON >= -140)]

        m1 = self.matches[(self.matches.SKS_PP_LON < -125) & (self.matches.SKS_PP_LON >= -140)]

        d2 = self.diffs[(self.diffs.SKS_PP_LON >= -125)]
        ns2 = self.null_split[(self.null_split.SKS_PP_LON >= -125)]
        m2 = self.matches[(self.matches.SKS_PP_LON >= -125) ]

        lat_diff1 = (d1.SKS_PP_LAT + d1.SKKS_PP_LAT)/2
        lat_match1 = (m1.SKS_PP_LAT + m1.SKKS_PP_LAT)/2
        lat_ns1 = (ns1.SKS_PP_LAT + ns1.SKKS_PP_LAT)/2

        lat_diff2 = (d2.SKS_PP_LAT + d2.SKKS_PP_LAT)/2
        lat_match2 = (m2.SKS_PP_LAT + m2.SKKS_PP_LAT)/2
        lat_ns2 = (ns2.SKS_PP_LAT + ns2.SKKS_PP_LAT)/2
        # Plot lat v lam2bar
        #For West half
        ax1.plot(d1.LAM2_BAR,lat_diff1,'k.',label= 'Discrepant')
        ax1.plot(ns1.LAM2_BAR,lat_ns1,'k.')
        ax1.plot(m1.LAM2_BAR,lat_match1,'kx',label='Matching')
        # For eastern half
        # ax3.plot(d2.LAM2_BAR,lat_diff2,'g.',label= 'Discrepant')
        # ax3.plot(ns2.LAM2_BAR,lat_ns2,'g.')
        # ax3.plot(m2.LAM2_BAR,lat_match2,'b.',label='Matching')
        #Plto lat v dSI
        ax2.plot(d1.D_SI_Pr,lat_diff1,'k.',label= 'Discrepant')
        ax2.plot(ns1.D_SI_Pr,lat_ns1,'k.',label='Discrepant')
        ax2.plot(m1.D_SI_Pr,lat_match1,'kx',label='Matching')
        # For Eastern half
        # ax4.plot(d2.D_SI_Pr,lat_diff2,'g.',label= 'Discrepant')
        # ax4.plot(ns2.D_SI_Pr,lat_ns2,'g.',label='Discrepant')
        # ax4.plot(m2.D_SI_Pr,lat_match2,'b.',label='Matching')
        ax1.set_xlim([0,0.1])
        ax2.set_xlim([0,1.2])
        ax1.set_ylim([36,54])
        ax1.set_ylabel(r'Latitude $(\degree)$')
        ax1.set_xlabel(r'$\bar{\lambda_2}$')
        ax2.set_xlabel(r'$\Delta SI$')
        ax2.set_ylabel(r'Latitude $(\degree)$')
        ax1.legend()
        plt.suptitle(r'Latitude v Discrepancy parameters for N-S trending SKS-SKKS pairs')
        plt.savefig('/Users/ja17375/DiscrePy/Figures/Lat_v_splitting_params.png',dpi = 400,transparent=True)
        plt.show()



#####################################################################################################
# Top level where script is invoked from command line
if __name__ == '__main__':
    print('Hello World, I am SDB_analysis.py and now its time to do my thing')

    if len(sys.argv) < 3:
        print('Excuse me, but you appear to have forgotten to provide enough arguements, 2 arguements (path and sdb_stem) are needed. Lets do this now.')
        path = input('Input the path to the Sheba/Results directory that contains the sdb files \n > ')
        sdb_stem = input('Input the stem of the SKS and SKKS sdb files you want to analse \n > ')


    elif len(sys.argv) == 3:
        path = sys.argv[1]
        sdb_stem = sys.argv[2]

    else:
        Warning('Too many arguements, script will not do anything else')

    print('{}/{}_SKS_SKKS.pairs'.format(path,sdb_stem))
    if os.path.isfile('{}/{}_SKS_SKKS.pairs'.format(path,sdb_stem)) is False:
        "Lets's make some pairs"
        B = Builder(p=path,RunDir=path.split('/')[-1],sdb_stem=sdb_stem,snr=1)
        B.run()
    else:
        'Pairs already exist'

    # # Now lets read the pairs
    # pair_file = '{}/{}_SKS_SKKS.pairs'.format(path,sdb_stem)
    # p = Pairs(pair_file)
    # p.match(pair_file[:-6])

    print('End')
