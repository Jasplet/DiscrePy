#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import sys
import os
import shlex
from subprocess import call
import matplotlib.pyplot as plt
import numpy as np
from stack import Stacker
from time import ctime
from glob import glob

class Builder:
    """
    Class to construct a .pairs Database (and Pierce Points), calc D_SI values and stack L2 surfaces.

    - sdb_stem: [str] sdb file stem for the data that you want to process
    - path: [str] path to Sheba/Results directory that contains the sdb's that you want
    - RunDir: [str] path to the Run Directory than Contains the .lamR surfaces to stack
    """
    def __init__(self,p,RunDir,sdb_stem,snr=2.0,syn=False):

        self.path = p
        self.path_stk = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/{}'.format(RunDir)
        self.sdb_stem = sdb_stem
        self.fpath =     '{}/{}_{:02d}_SKS_SKKS.pairs'.format(self.path,self.sdb_stem,int(snr))
        self.lam2 = [ ]
        self.snr = snr
        self.syn=syn
        #if kwargs are none:\

    def run(self):
        '''Function to build the pairs file'''
        # First match the SKS and SKKS
        start = ctime()
        print('Making Pairs')
        self.make_pairs()
        # Apply a quick Signal to Noise test to get rid of the rreally bad data
        self.P = self.snr_check() # Overwrite self.P weith only accepted events
        self.write_out(self.P,name='{}_all.pairs'.format(self.sdb_stem))
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
            self.match()
            print(r'Apply $\bar{\lambda_2}$ & $\Delta SI$ test')
            self.match_l2()
            # print('{} pairs'.format(len(self.P)))

        # And save the result
        self.write_out(self.P,'{}_{:02d}.pairs'.format(self.sdb_stem,int(self.snr)))
        end = ctime()
        print('END. start {}, end {}'.format(start,end))

    def gen_pp(self):
        ''' Fucntion to test for whether the .pp file exists and if not call TauP to generate it and the corresponding mspp files '''

        print('Looking for pp file {}.pp'.format(self.fpath.split('.')[0]))
        if os.path.isfile('{}.pp'.format(self.fpath.split('.')[0])):
            #print(pf[:-6])
            self.pp = pd.read_csv('{}.pp'.format(self.fpath.split('.')[0]),delim_whitespace=True)
        else:
            print('Pierce Points file {}.pp doesnt not exist, calling pierce.sh'.format(self.fpath.split('.')[0]))
            p = call(shlex.split('/Users/ja17375/Shear_Wave_Splitting/Sheba/Programs/pierce.sh {}'.format(self.fpath)))
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
              'FAST_y':'FAST_SKKS', 'DFAST_y': 'DFAST_SKKS','TLAG_y':'TLAG_SKKS','DTLAG_y':'DTLAG_SKKS','SPOL_y':'SPOL_SKKS','DSPOL_y':'DSPOL_SKKS',
              'WBEG_y':'WBEG_SKKS','WEND_y':'WEND_SKKS','EIGORIG_y':'EIGORIG_SKKS','EIGCORR_y':'EIGCORR_SKKS','Q_y':'Q_SKKS','SNR_y':'SNR_SKKS','NDF_y':'NDF_SKKS'}
        # The dictionary relabels the other columns in the join so that we can more easily pick apart the SKS and SKKS results
        SKS_SKKS_pair.rename(relabel,axis='columns',inplace=True)
        # Sort the Pairs dataframe so the pairs are in chronological order (by origin time (DATE only))
        self.P = SKS_SKKS_pair.sort_values(by=['DATE'],ascending=True)

    def write_out(self,df,name):
        df.to_csv('{}/{}'.format(self.path,name),sep=' ',index=False)

    def add_DSI(self):
        '''Calculate the difference in Splitting Intensity for each pair and add it to dataframe'''
        si_sks = self.P.INTENS_x
        si_skks = self.P.INTENS_y
        d_si = np.abs(si_sks-si_skks)
        self.P['D_SI'] = d_si
        #Delete SI cols as we dont need them any more ?
        del self.P['INTENS_x']
        del self.P['INTENS_y']

    def pair_stack(self,mode='man'):
        ''' Runs Stacker for all the desired pairs (a .pairs file)'''

        ext ='lamR'
        print('Stacking')
        rd = self.path_stk.split('/')[-1]
        out = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}/Stacks'.format(rd)
        if os.path.isdir(out) is False:
            print('{} does not exist, creating'.format(out))
            os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/{}/Stacks'.format(rd))

        for i,f in enumerate(self.P.DATE.values):
            # print('It {}, time is {} '.format(i,str(datetime.now())))
            # First get the right DATE,TIME and STATION
            date,time,stat = self.P.DATE[i], self.P.TIME[i], self.P.STAT[i]
            fstem = '{}_{}_{}'.format(stat,date,time)

            lam2_stem = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path_stk,stat,fstem,ext))
            # print(lam2_stem)
            print('{}/{}/SKS/{}??_SKS.lam2'.format(self.path_stk,stat,fstem))
            if len(lam2_stem) is not 0:
                # I.e if glob has managed to find the sks lam2 surface file
                sks_lam2 = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path_stk,stat,fstem,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.{}'.format(self.path_stk,stat,fstem,ext))[0]
                Stk = Stacker(sks_lam2,skks_lam2,out)
                if mode == 'man':
                    self.lam2.append(Stk.lam2)
                elif mode == 'sheba':
                    self.lam2.append(Stk.sol[-1])
            else:
                fstem2 = '{}_{}'.format(stat,date)
                print('fstem2')
                print('{}/{}/SKS/{}_*_SKS.{}'.format(self.path_stk,stat,fstem2,ext))
                sks_lam2 = glob('{}/{}/SKS/{}_*_SKS.{}'.format(self.path_stk,stat,fstem2,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}_*_SKKS.{}'.format(self.path_stk,stat,fstem2,ext))[0]
                # Now for a sanity check
                if (len(sks_lam2) is not 0) or (len(skks_lam2) is not 0):
                    Stk = Stacker(sks_lam2,skks_lam2,out)
                    if mode == 'man':
                        self.lam2.append(Stk.lam2)
                    elif mode == 'sheba':
                        self.lam2.append(Stk.sol[-1])
                else:
                    #print('lam2 surfaces cannot be found, skipping')
                    pass

    def add_lam2(self):
        '''
        Stack the associated raw lambda 2 surfaces (as output by sheba) for each SKS SKKS pair and find min value
        '''

        self.pair_stack()
        l2df = {'LAM2' : self.lam2}
        ldf = pd.DataFrame(l2df)
        self.P['LAM2'] = ldf

    def add_pp(self):
        '''Adds piercepoints to .pairs file'''
        if len(self.pp) == len(self.P):
            print(len(self.pp) ,len(self.P))
            self.P['SKS_PP_LAT'] = self.pp.lat_SKS
            self.P['SKS_PP_LON'] = self.pp.lon_SKS
            self.P['SKKS_PP_LAT'] = self.pp.lat_SKKS
            self.P['SKKS_PP_LON'] = self.pp.lon_SKKS
        else:
            print('Dimension mismatch, not adding piercepoints')
    def match(self,sigma=2):
        """
        Funntion to see if the SKS and SKKS splititng measurements for a pair of measurements match within error
        file - the filename you want for the output files
        Default error for this kind of anlysis is 2-sigma. Sheba returns 1 sigma so the DFAST and DTLAG need to be scaled appropriatly.
        """
        fstem = self.path.split('.')[0]# Split('.')[0] takes off the extension
        # First Add piercepoints to pairas files
        self.add_pp()
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

    def match_l2(self,t_l2 = 0.01, t_dSI = 0.15):
        """
        Function to test if a SK(K)S pair are matching or discrepant. This test is based on lb2 and dSI (see synthetics_stakcs notes)
        This function assumes that lam2 and dSI have already been calculated and added to the pairs file. It also assumes that you
        want to keep the same file names but just add _match_l2 and _diff_l2self.
        The lam2 threshold used is 0.01 and the dSI threshold is 0.15
        """

        fstem = self.path.split('.')[0]# Split('.')[0] takes off the extension
        # Now apply the test to find the discrepant pairs, by definition the remainder must by the matches
        null_pairs  = self.P[((self.P.Q_SKS < -0.5) & (self.P.Q_SKKS < -0.5))] # Pairs where both phases are nulls (according to Q), auto classify as matching
        null_split_pair = self.P[(((self.P.Q_SKS < -0.5) & (self.P.Q_SKKS > 0.5)) | ((self.P.Q_SKS > 0.5) & (self.P.Q_SKKS < -0.5)))] # Test for pairs with 1 null 1 split, discrepant by definition
        splits = self.P[((self.P.Q_SKS > 0.5) & (self.P.Q_SKKS > 0.5 ))] # Test for pairs whjere both phases are split
        s_diff= splits[((self.P.LAM2 >= t_l2) &  (self.P.D_SI > t_dSI))] # Apply tests for discrepant splitting
        s_match = splits[((self.P.LAM2 < t_l2) & (self.P.D_SI < t_dSI))]

        # Now combined matching and discrepant pairs together
        diff = null_split_pair.append(s_diff) # Combine the null-split pairs and pairs of discrepant splitting
        match = null_pairs.append(s_match)
        uID_int = self.P.drop(diff.index)
        uID     = uID_int.drop(match.index)
        match.to_csv('{}/{}_{:02d}_matches_l2.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        diff.to_csv('{}/{}_{:02d}_diffs_l2.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        uID.to_csv('{}/{}_{:02d}_uDI_l2.pairs'.format(self.path,self.sdb_stem,int(self.snr)),index=False,sep=' ')
        # Open up mspp files
        print('Writing to {}/{}_{:02d}'.format(self.path,self.sdb_stem,int(self.snr)))
        mspp_match = open('{}/{}_{:02d}_matches_l2.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_diff = open('{}/{}_{:02d}_diffs_l2.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')
        mspp_uID = open('{}/{}_{:02d}_diffs_l2.mspp'.format(self.path,self.sdb_stem,int(self.snr)),'w+')

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



    def snr_check(self):
        ''' Run a quick and dirty check on Signal to Noise Ratio. If is is less than the threshold events are rejected'''
        t = self.snr
        print(t)
        self.accepted_i = [ ]
        self.d = [ ]
        print('There are {} pairs pre-SNR < 2 test'.format(len(self.P)))
        for i,row in self.P.iterrows():
            if row.SNR_SKS <= t or row.SNR_SKKS <= t:
                #Test to see if Signal-to-Noise is too high
                print('SNR for SKS or SKKS less than {:02}, auto-reject'.format(t))
                self.d.append(i)
            else:
                self.accepted_i.append(i)
                # print('Event accepted')

        print('{} accepted, {} rejected'.format(len(self.accepted_i),len(self.d)))
        accepted_pairs = self.P.drop(self.d)
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
                # Read in standard match,diff files
                pmatch = '{}_matches.pairs'.format(fname.split('.')[0])
                pdiff ='{}_diffs.pairs'.format(fname.split('.')[0])
                self.matches = pd.read_csv(pmatch,delim_whitespace=True,converters=date_time_convert)
                self.diffs = pd.read_csv(pdiff,delim_whitespace=True,converters=date_time_convert)
                # Now read in new match,diffs made using lam2/dSI test
                pmatch = '{}_matches_l2.pairs'.format(fname.split('.')[0])
                pdiff ='{}_diffs_l2.pairs'.format(fname.split('.')[0])
                self.matches_l2 = pd.read_csv(pmatch,delim_whitespace=True,converters=date_time_convert)
                self.diffs_l2 = pd.read_csv(pdiff,delim_whitespace=True,converters=date_time_convert)
            except FileNotFoundError:
                print('Match or Diff file not found. Are you using synhetics maybe??')


        elif file is False:
            print('Expectinf df input')
            self.df = df

    def plot_dist_v_discrep(self):

        fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2, figsize = [12,6])

        ax1.plot(self.matches.DIST,self.matches.LAM2,marker='.',color='blue',label='Matching')
        ax1.plot(self.diffs.DIST,self.diffs.LAM2,marker='.',color='darkorange',label='Discrepant')
        ax1.set_xlim([105, 140])
        ax1.set_ylim([0, np.around(self.df.LAM2.max(),decimals=1)])
        ax1.set_xlabel('Epicentral Distance (Deg)')
        ax1.set_ylabel(r'$\bar{\lambda _2}$ values')
        ax1.set_title('Match/Diff according to  2 sigma test')
        ax1.legend()

        ax2.plot(self.matches.DIST,self.matches.D_SI,marker='.',color='blue',label='Matching')
        ax2.plot(self.diffs.DIST,self.diffs.D_SI,marker='.',color='darkorange',label='Discrepant')
        ax2.set_xlim([105,140])
        ax2.set_ylim([0,np.around(self.df.D_SI.max(),decimals=1)])
        ax2.set_xlabel('Epicentral Distrance (Deg)')
        ax2.set_ylabel(R'$\Delta SI $')

        ax3.plot(self.matches_l2.DIST,self.matches_l2.LAM2,marker='.',color='blue',label='Matching')
        ax3.plot(self.diffs_l2.DIST,self.diffs_l2.LAM2,marker='.',color='darkorange',label='Discrepant')
        ax3.set_xlim([105, 140])
        ax3.set_ylim([0, np.around(self.df.LAM2.max(),decimals=1)])
        ax3.set_xlabel('Epicentral Distance (Deg)')
        ax3.set_ylabel(r'$\bar{\lambda _2}$ values')
        ax3.set_title(r'Match/Diff according to $\bar{lambda _2} / \Delta SI $ test')
        ax3.legend()

        ax4.plot(self.matches_l2.DIST,self.matches_l2.D_SI,marker='.',color='blue',label='Matching')
        ax4.plot(self.diffs_l2.DIST,self.diffs_l2.D_SI,marker='.',color='darkorange',label='Discrepant')
        ax4.set_xlim([105,140])
        ax4.set_ylim([0,np.around(self.df.D_SI.max(),decimals=1)])
        ax4.set_xlabel('Epicentral Distrance (Deg)')
        ax4.set_ylabel(R'$\Delta SI $')

        plt.show()

    def lam2_v_SI(self,figname=None,save=False):
        '''Plot a scatter plot of lambda2 values verus splitting INtensity difference.For both methods of categorising matching and discrepanct results'''

        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(16,8))
        # Isloate clear Splits
        ms_sig= self.matches[(self.matches.Q_SKS > 0.7) & (self.matches.Q_SKKS > 0.7)]
        ds_sig = self.diffs[(self.diffs.Q_SKS > 0.7) & (self.diffs.Q_SKKS > 0.7)]
        ms_lam = self.matches_l2[(self.matches_l2.Q_SKS > 0.7) & (self.matches_l2.Q_SKKS >0.7)]
        ds_lam = self.diffs_l2[(self.diffs_l2.Q_SKS > 0.7) & (self.diffs_l2.Q_SKKS >0.7)]
        # Isolate clear nulls
        mn_sig = self.matches[(self.matches.Q_SKS < -0.7) & (self.matches.Q_SKKS < -0.7)]
        dn_sig = self.diffs[(self.diffs.Q_SKS < -0.7) & (self.diffs.Q_SKKS < -0.7)]
        mn_lam = self.matches_l2[(self.matches_l2.Q_SKS < -0.7) & (self.matches_l2.Q_SKKS < -0.7)]
        dn_lam = self.diffs_l2[(self.diffs_l2.Q_SKS < -0.7) & (self.diffs_l2.Q_SKKS < -0.7)]
        # Plot clear splits - matches and diff according to 2 sigma
        ax1.plot(ms_sig.LAM2,ms_sig.D_SI,color='blue',marker='.',ls='None',label="'Matching' splits")
        ax1.plot(ds_sig.LAM2,ds_sig.D_SI,color='darkorange',marker='.',ls='None',label="'Discrepant' splits")
        # # Plot clear nulls - matches and diff according to 2 sigma
        # ax1.plot(mn_sig.LAM2,mn_sig.D_SI,color='blue',marker='x',ls='None',label="'Matching' nulls")
        # ax1.plot(dn_sig.LAM2,dn_sig.D_SI,color='darkorange',marker='x',ls='None',label="'Discrepant' nulls")
        # Plot Deng threshold
        # ax1.plot([0, 1.0],[0.4,0.4],ls='dashed',color='black')
        ax1.set_xlabel(r'$\bar{\lambda _2}$')
        ax1.set_ylabel(r'$\Delta$ SI')
        ax1.set_title(r'Matching/Discrepant pairs according to $2 \sigma$.')
        ax1.legend()
        # ax1.set_ylim([0,max([np.around(ms_sig.D_SI.max(),decimals=1),np.around(mn_sig.D_SI.max(),decimals=1)])])
        # ax1.set_xlim([0,max([np.around(ms_sig.LAM2.max(),decimals=1),np.around(mn_sig.LAM2.max(),decimals=1)])])
        #Lam2 / dSI test
        # Plot clear splits - matches and diff according to lam2/dSI
        ax2.plot(ms_lam.LAM2,ms_lam.D_SI,color='blue',marker='.',ls='None',label="'Matching' splits")
        ax2.plot(ds_lam.LAM2,ds_lam.D_SI,color='darkorange',marker='.',ls='None',label="'Discrepant' splits")
        # Plot clear nulls - matches and diff according to lam2/dSI
        # ax2.plot(mn_lam.LAM2,mn_lam.D_SI,color='blue',marker='x',ls='None',label="'Matching' nulls")
        # ax2.plot(dn_lam.LAM2,dn_lam.D_SI,color='darkorange',marker='x',ls='None',label="'Discrepant' nulls")
        # ax2.plot([0, 1.0],[0.4,0.4],ls='dashed',color='black')
        ax2.set_xlabel(r'$\bar{\lambda _2}$')
        ax2.set_ylabel(r'$\Delta$ SI')
        ax2.legend()
        ax2.set_title(r'Matching/Discrepant pairs according to $\bar{\lambda _2}$ & $\Delta SI$.')
        # ax2.set_ylim([0,max([np.around(ms_lam.D_SI.max(),decimals=1),np.around(mn_lam.D_SI.max(),decimals=1)])])
        # ax2.set_xlim([0,max([np.around(ms_lam.LAM2.max(),decimals=1),np.around(mn_lam.LAM2.max(),decimals=1)])])
        if save is True:
            plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/{}.eps'.format(figname),format='eps',dpi=1000)

        plt.show()

    def l2_dSI_hist(self):

        fig,(ax1,ax2) = plt.subplots(1,2,figsize=(14,6))
        bins_l2 = np.linspace(0,0.5,10)
        bins_dsi = np.linspace(0,3.5,10)
        m_splits = self.matches_l2[(self.matches_l2.Q_SKS > 0.5) & (self.matches_l2.Q_SKKS > 0.5)]
        d_splits = self.diffs_l2[(self.diffs_l2.Q_SKS > 0.5) & (self.diffs_l2.Q_SKKS > 0.5)]

        ax1.hist([m_splits.LAM2,d_splits.LAM2],bins_l2, histtype='bar', stacked=True,label=["'Matching'","'Discrepent'"])
        ax1.set_xlabel(r'$\bar {\lambda _2}$ values')
        ax1.set_ylabel('Count')
        ax1.legend()
        ax1.set_xlim([0, np.around(self.df.LAM2.max(),decimals=1)])
        ax2.hist([m_splits.D_SI,d_splits.D_SI],bins_dsi, histtype='bar', stacked=True,label=["'Matching'","'Discrepent'"])
        ax2.set_xlabel(r'$\Delta SI$ values')
        ax2.set_ylabel('Count')
        ax2.legend()
        ax2.set_xlim([0,np.around(self.df.D_SI.max(),decimals=1)])
        # ax2.set_ylim([0,4.0])
        plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/hist_and_dVs.eps',format='eps',dpi=1000)
        plt.show()

    def plot_dist_v_split(self,save=False):

        # Plot colored by Lambda 2
        fig1,(ax1,ax2) = plt.subplots(1,2, figsize = [12,6])

        deltaF = 90 -abs(abs(self.df.FAST_SKS - self.df.FAST_SKKS) - 90 ) # This should allow for differences to vary between 0 and 90 whilst also dealing with the wrap around properly
        ddt = abs(self.df.TLAG_SKS) - abs(self.df.TLAG_SKKS)

        C1 = ax1.scatter(self.df.DIST,deltaF,c=self.df.LAM2,marker='.',label=r'$| \phi_{SKS} - \phi_{SKKS} |')

        ax1.set_xlim([105, 140])
        ax1.set_ylim([0,90])
        ax1.set_xlabel('Epicentral Distance (Deg)')
        ax1.set_ylabel(r'$ |\phi_{SKS} - \phi_{SKKS}|$')
        ax1.set_title(r'Difference in $\phi$, coloured by $\Lambda$')

        C2 = ax2.scatter(self.df.DIST,ddt,c=self.df.LAM2,marker='.',label='delta deltaT')
        #ax2.plot(df.diffs.DIST,df.diffs.D_SI,'r.',label='Discrepant')
        ax2.set_xlim([105,140])
        ax2.set_ylim([0,4])
        ax2.set_xlabel('Epicentral Distrance (Deg)')
        ax2.set_ylabel(r'$| \delta t_{SKS} - \delta t_{SKKS} | $')
        ax2.set_title(r'Difference in $(\delta_t)$, coloured by $\Lambda$')
        cbar1 = fig1.colorbar(C1,use_gridspec=True)
        cbar1.set_label(r'$\Lambda$',rotation=270)

        fig2,(ax3) = plt.subplots(figsize= [6,6])
        C3 = ax3.scatter(ddt,(deltaF),c=self.df.LAM2,marker='.')
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
            fig1.savefig('/Users/ja17375/Shear_Wave_Splitting/Sheba/Figures/diff_v_dist_Lambda.eps',format='eps',transparent=True,dpi=400)
            fig2.savefig('/Users/ja17375/Shear_Wave_Splitting/Sheba/Figures/delfast_v_deltlag_Lambda.eps',format='eps',transparent=True,dpi=400)
            fig3.savefig('/Users/ja17375/Shear_Wave_Splitting/Sheba/Figures/diff_v_dist_DSI.eps',format='eps',transparent=True,dpi=400)
            fig4.savefig('/Users/ja17375/Shear_Wave_Splitting/Sheba/Figures/delfast_v_deltlag_DSI.eps',format='eps',transparent=True,dpi=400)

            plt.close('all')
        else:
            plt.show()

    def Q_v_l2_dSI(self):
        '''
        Make plots of Q factor (Wuestefeld et al) sgainst lambda 2 and Delta SI
        '''
        fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(14,6),sharex=True)#,sharey=True)

        C1 = ax1.scatter(self.df.Q_SKS,self.df.LAM2,c=self.df.SNR_SKS,marker='.',vmin=2,vmax=20)
        # ax1.plot(self.matches_l2.Q_SKS,self.matches_l2.LAM2,marker='.',color='blue',ls='None')
        #ax1.plot(self.diffs_l2.Q_SKS,self.diffs_l2.LAM2,marker='.',color='darkorange',ls='None')
        ax1.set_xlabel('SKS Q factor')
        ax1.set_ylabel(r'$\bar{\lambda _2}$')
        ax1.set_xlim([-1,1])
        ax1.set_ylim([0,np.around(self.df.LAM2.max(),decimals=1)])
        plt.colorbar(C1,ax=ax1)

        ax2.plot(self.matches_l2.Q_SKKS,self.matches_l2.LAM2,marker='.',color='blue',ls='None')
        #ax2.plot(self.diffs_l2.Q_SKKS,self.diffs_l2.LAM2,marker='.',color='darkorange',ls='None')
        ax2.set_xlabel('SKKS Q factor')
        ax2.set_xlim([-1,1])
        ax2.set_ylim([0,np.around(self.df.LAM2.max(),decimals=1)])

        ax3.plot(self.matches_l2.Q_SKS,self.matches_l2.D_SI,marker='.',color='blue',ls='None')
        #ax3.plot(self.diffs_l2.Q_SKS,self.diffs_l2.D_SI,marker='.',color='darkorange',ls='None')
        ax3.set_xlabel('SKS Q factor')
        ax3.set_ylabel(r'$\Delta SI$')
        ax3.set_xlim([-1,1])
        ax3.set_ylim([0,np.around(self.df.D_SI.max(),decimals=1)])

        ax4.plot(self.matches_l2.Q_SKKS,self.matches_l2.D_SI,marker='.',color='blue',ls='None')
        #ax4.plot(self.diffs_l2.Q_SKKS,self.diffs_l2.D_SI,marker='.',color='darkorange',ls='None')
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
        ax1.hist(self.df.SNR_SKS,bins=np.arange(0,50,5),histtype='bar')
        ax1.set_xlabel('S/N ratio (SKS)')
        ax1.set_ylabel('Count')

        ax2.hist(self.df.SNR_SKKS,bins=np.arange(0,50,5),histtype='bar')
        ax2.set_xlabel('S/N ratio (SKKS)')
        ax2.set_ylabel('Count')
        plt.show()


    def plot_SNR_v_l2(self):
        '''
        Make plots of SNR v lambda2. In the style of phi_i v SNR from Restivo and Helffrich (2006)
        '''
        fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,sharey='row',figsize=(8,8))
        #Plot SNR for SKS
        ax1.plot(self.df.SNR_SKS,self.df.LAM2,'k.')
        ax1.set_ylabel(r'$\bar{\lambda_2}$')
        ax1.set_xlabel('S/N ratio')
        ax1.set_title('d fast SKS determination dependance on S/N')
        #Plot SNR for SKKS
        ax2.plot(self.df.SNR_SKKS,self.df.LAM2,'k.')
        ax2.set_ylabel(r'$\bar{\lambda_2}$')
        ax2.set_xlabel('S/N ratio')
        ax2.set_title('d fast SKKS determination dependance on S/N')
        plt.show()

    def package(self,outdir,mypath='/Users/ja17375/Shear_Wave_Splitting/Sheba/SAC'):
        '''
        Function to package the SAC files for the observations made into a new directory (for sharing data and results)
        '''
        if os.path.isdir('{}{}'.format(mypath,outdir)) is False:
            #If outdir doesnt exist then Make it ad underlying structure (an additional SAC directory )
            os.makedirs('{}{}/SAC'.format(mypath,outdir))

        def mk_sacfiles(path,stat,date,time):
            '''Function to generate sacfile names'''
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
        self.p_sorted = self.df.sort_values(by='LAM2',ascending=True)
        self.p_sorted.reset_index(drop=True)

        # Find indecies of events we want to plot
        if surfs_to_plot is None:
            self.surfs= np.round(np.arange(0,len(self.p_sorted),round((len(self.p_sorted)/nplots))))
        else:
            self.surfs = list(set([i if i < len(self.p_sorted) else (len(self.p_sorted)-1) for i in surfs_to_plot])) # This makes sure that indicies are always within the range of available surfaces (stops errors for occuring)
            #self.surfs = self.surfs.sort() # Sorts list in ascending order, has to be done speratly as sort acts of list and returns nothing
        print(self.surfs)
        if save is True:
            dir = input('Enter Directory you want to save stacked surfaces to > ')

            if os.path.isdir('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}'.format(dir)) is False:
                make_d = input('Directory {} does not exist, do you want to create it? (y/n)'.format(dir))
                if make_d =='y':
                    print('Ok, creating directory {}'.format(dir))
                    os.mkdir('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}'.format(dir))
                else:
                    print('Exiting....')
                    sys.exit()

        for s in self.surfs:
            print(s)
            stat,date,time = self.p_sorted.STAT.values[s],self.p_sorted.DATE.values[s], self.p_sorted.TIME.values[s]
            # Note that dtlag and dfast are multiplied through by sigma HERE !
            if self.syn is True:
                self.lam2_surface(f1=self.syn1[s],f2=self.syn2[s])
                fast_sks,dfast_sks,lag_sks,dlag_sks = self.df.FAST_SKS.values[s],(sigma*self.df.DFAST_SKS.values[s]),self.df.TLAG_SKS.values[s],(sigma*self.df.DTLAG_SKS.values[s])
                fast_skks,dfast_skks,lag_skks,dlag_skks = self.df.FAST_SKKS.values[s],(sigma*self.df.DFAST_SKKS.values[s]),self.df.TLAG_SKKS.values[s],(sigma*self.df.DTLAG_SKKS.values[s])
                lam2 = self.df.LAM2.values[s]
            else:
                fast_sks,dfast_sks,lag_sks,dlag_sks = self.p_sorted.FAST_SKS.values[s],(sigma*self.p_sorted.DFAST_SKS.values[s]),self.p_sorted.TLAG_SKS.values[s],(sigma*self.p_sorted.DTLAG_SKS.values[s])
                fast_skks,dfast_skks,lag_skks,dlag_skks = self.p_sorted.FAST_SKKS.values[s],(sigma*self.p_sorted.DFAST_SKKS.values[s]),self.p_sorted.TLAG_SKKS.values[s],(sigma*self.p_sorted.DTLAG_SKKS.values[s])
                lam2 = self.p_sorted.LAM2.values[s]
                l_path = '{}/{}_{}_{}'.format(self.spath,stat,date,time) #Path to lambda 2 surfaces for SKS and SKKS
                print(l_path)
                self.lam2_surface(l_path)

            print('Stat {}, Evt Time {}-{} LAM2 = {}'.format(stat,date,time,lam2))
            # fig = plt.figure(figsize=(12,12))
            fig, (ax0,ax1,ax2) = plt.subplots(1,3,figsize=(24,8),sharey=True)
            fig.patch.set_facecolor('None')
            if self.syn is True:
                plt.suptitle(r'Syn Stack E1: {} E2: {} $\lambda _2$ value = {:4.3f}'.format(self.syn1[s],self.syn2[s],self.df.LAM2.values[s]),fontsize=28)
            else:
                plt.suptitle(r'Event {}_{}_{}. Stacked $\lambda _2$ value = {:4.3f}'.format(stat,date,time,self.p_sorted.LAM2.values[s]),fontsize=28)

            # gs = gridspec.GridSpec(3,2)
            # ax0 = plt.subplot(gs[0,0])
            ax0.set_title(r'SKS $\lambda _2$ surface',fontsize=24)
            C0 = ax0.contourf(self.T,self.F,(self.sks_lam2),[0.005,0.01,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5],cmap='inferno_r',extend='both')
            ax0.contour(C0,colors='k')
            # ax0.clabel(C0,C0.levels,inline=True,fmt ='%2.3f')
            #Plot SKS Solution
            ax0.plot(lag_sks,fast_sks,'b.',label='SKS Solution')
            print('Lag sks {}. Fast SKS {}.'.format(lag_sks,fast_sks))
            ax0.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax0.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            ax0.set_ylabel(r'Fast,$\phi$, (deg)')
            ax0.set_xlabel(r'Lag ,$\delta$ t, (sec)')
            #Plot SKKS Solution
            ax0.plot(lag_skks,fast_skks,'r.',label='SKKS Solution')
            ax0.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'r-')
            ax0.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'r-')
            ax0.set_ylim([-90,90])
            ax0.set_xlim([0,4])
            ax0.set_yticks([-90,-60,-30,0,30,60,90])
            # ax0.contourf(self.sks_lam2,cmap='inferno_r')
            # ax1 = plt.subplot(gs[0,1])
            C1 = ax1.contourf(self.T,self.F,self.skks_lam2,[0,0.005,0.01,0.025,0.05,0.075,0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.25,1.5],cmap='inferno_r',extend='both')
            C3 = ax1.contour(C1,colors='k')
            # ax1.clabel(C1,C1.levels,inline=True,fmt ='%2.3f')
            # ax1.contourf(self.skks_lam2,cmap='magma')
            #Plot SKS Solution
            ax1.plot(lag_sks,fast_sks,'b.',label='SKS Solution')
            ax1.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax1.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            #Plot SKKS Solution
            ax1.plot(lag_skks,fast_skks,'r.',label='SKKS Solution')
            ax1.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'r-')
            ax1.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'r-')
            ax1.set_xlabel(r'Lag ,$\delta$ t, (sec)')
            ax1.set_ylim([-90,90])
            ax1.set_xlim([0,4])
            ax1.set_yticks([-90,-60,-30,0,30,60,90])
            ax1.set_title(r'SKKS $\lambda _2$ surface',fontsize=24)
            # ax2 = plt.subplot(gs[1:,:])
            if self.syn is True:
                self.show_stacks(ax2,evt=self.synstk[s],path='/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/Stacks')
            else:
                self.show_stacks(ax2,l_path.split('/')[-1])
            # print('STK_FAST: {} +/- {}'.format(self.df_fast,self.df_dfast))
            # Modify stk_dlag and stk_dfast by sigma
            ##########################
            # self.stk_dlag = self.stk_dlag*sigma
            # self.stk_dfast = self.stk_dfast*sigma
            ############################
            ax2.set_ylim([-90,90])
            ax2.set_xlim([0,4])
            ax2.set_yticks([-90,-60,-30,0,30,60,90])
            #Plto SKS solution
            ax2.plot(lag_sks,fast_sks,'b.',label='SKS Solution')
            ax2.plot([lag_sks-dlag_sks,lag_sks+dlag_sks],[fast_sks,fast_sks],'b-')
            ax2.plot([lag_sks,lag_sks],[fast_sks-dfast_sks,fast_sks+dfast_sks],'b-')
            #Plot SKKS solution
            ax2.plot(lag_skks,fast_skks,'r.',label='SKKS Solution')
            ax2.plot([lag_skks-dlag_skks,lag_skks+dlag_skks],[fast_skks,fast_skks],'r-')
            ax2.plot([lag_skks,lag_skks],[fast_skks-dfast_skks,fast_skks+dfast_skks],'r-')
            # Plot Stacked solution on SKS
            ax0.plot(self.stk_lag,self.stk_fast,'g.',label='Stacked Solution')
            # ax0.plot([self.stk_lag- self.stk_dlag, self.stk_lag + self.stk_dlag],[ self.stk_fast, self.stk_fast],'g-')
            # ax0.plot([self.stk_lag, self.stk_lag],[ self.stk_fast - self.stk_dfast ,self.stk_fast +self.stk_dfast],'g-')
            # Plot Stacked Solution on SKKS surface
            ax1.plot(self.stk_lag,self.stk_fast,'g.',label='Stacked Solution')
            # ax1.plot([self.stk_lag-self.stk_dlag,self.stk_lag+self.stk_dlag],[self.stk_fast,self.stk_fast],'g-')
            # ax1.plot([self.stk_lag,self.stk_lag],[self.stk_fast-self.stk_dfast,self.stk_fast+self.stk_dfast],'g-')
            ## Add a legend (on ax0)
            ax0.legend(bbox_to_anchor=(0,1),loc='upper left')


            # cb = fig.colorbar(C1)
            # cb.add_lines(C3)

            ax2.set_title('Stacked SKS SKKS surface',fontsize=24)
            if save is True:
                # dir = input('Enter Directory you want to save stacked surfaces to > ')
                plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/Stacked_Surfaces/{}/LAM2_{:4.4f}_STAT_{}.eps'.format(dir,lam2,stat),format='eps',dpi=800)
                plt.close()
            elif save is False:
                plt.show()

    def lam2_surface(self,fstem=None,f1=None,f2=None):
        ''' Function to read  SKS and SKKS .lam2 surface files from sheba
        If syn if False (i.e real data is being used.) Then fstem in needed
        IF syn is True then f1 , f2 are needed
        '''
        # print(fstem)
        if self.syn == False:
            t_sks = '{}??_SKS.lamR'.format(fstem)
            t_skks = '{}??_SKKS.lamR'.format(fstem)

            sks =glob(t_sks)
            skks = glob(t_skks)
            print(sks)
            if len(sks) == 0:
                print('{}/SKS/{}*_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),stem))
                stem = '_'.join(fstem.split('/')[-1].split('_')[0:-1]) # aka cut off time part of file extension
                # print('{}/SKS/{}*_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),fstem.split('/')[-1]))
                sks = glob('{}/SKS/{}*_SKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),stem))
                skks = glob('{}/SKKS/{}*_SKKS.lamR'.format('/'.join(fstem.split('/')[0:-1]),stem))
                print(sks)

            self.sks_lam2 = np.loadtxt(sks[0])#,skiprows=4) # skip rows not needed for .lamR files
            self.skks_lam2 = np.loadtxt(skks[0])#,skiprows=4)

        elif self.syn == True:
            self.sks_lam2 = np.loadtxt('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/Stacks/{}'.format(f1))
            self.skks_lam2 = np.loadtxt('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SYNTH/Stacks/{}'.format(f2))

        nfast,nlag = self.sks_lam2.shape ;
        lag_max = 4.
        [self.T,self.F] = np.meshgrid(np.linspace(0,lag_max,num=nlag),np.arange(-90,91,1)) ;

    def show_stacks(self,ax,evt,path='/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Nulls'):
        '''Function to find and plot desired surface stacks based on the LAMDA2 value '''
        ### Plot Min Lamnda 2
        if self.syn is True:
            stk = np.loadtxt('{}/{}'.format(path,evt))
        else:
            print('{}/Stacks/{}??.lamSTK'.format(path,evt))
            file =glob('{}/Stacks/{}??.lamSTK'.format(path,evt))
            if len(file) == 0:
                file = glob('{}/Stacks/{}*.lamSTK'.format(path,'_'.join(evt.split('_')[0:-1])))

            stk = np.loadtxt(file[0])

        nfast,nlag = stk.shape ;
        lag_step = 0.025
        lag_max = (nlag) * lag_step;
        [T,F] = np.meshgrid(np.arange(0,lag_max,lag_step),np.arange(-90,91,1)) ;
        jf,jt  = np.unravel_index(stk.argmin(),stk.shape)
        print(stk.argmin())
        print(jf,jt)

        fast = np.arange(-90,91,1)[jf]
        lag = np.arange(0,lag_max,lag_step)[jt]
        C = ax.contourf(T,F,stk,[0.01,0.025,0.05,0.075,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,0.9,1.0,1.5],cmap='inferno_r',extend='both')
        C2 = ax.contour(C,colors='k')
        ax.set_ylabel(r'Fast,$\phi$, (deg)')
        ax.set_xlabel(r'Lag ,$\delta$ t, (sec)')
        # ax.plot([lag-dlag,lag+dlag],[fast,fast],'g-')
        # ax.plot([lag,lag],[fast-dfast,fast+dfast],'g-')
        ax.plot(lag,fast,'g.')
        # ax.clabel(C,C.levels,inline=True,fmt ='%4.3f')

        # self.cbar = plt.colorbar(C)
        # self.cbar.add_lines(C2)
        # ax.set_title(r'Event {}. $\lambda$ 2 value = {}'.format(evt,lam2))
        # Add fast, lag as attributes so we can plot them elsewhere
        sol = stk.min()

        print('Lam2 {}, fast {} lag {}'.format(sol,fast,lag))
        self.stk_fast = fast
        self.stk_lag = lag

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
        "Lets's find some pairs"
        make_pairs(path,sdb_stem)
    else:
        'Pairs already exist'

    # Now lets read the pairs
    pair_file = '{}/{}_SKS_SKKS.pairs'.format(path,sdb_stem)
    p = Pairs(pair_file)
    p.match(pair_file[:-6])

    print('End')
