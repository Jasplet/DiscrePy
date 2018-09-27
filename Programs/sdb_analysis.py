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
    def __init__(self,p,RunDir,sdb_stem):

        self.path = p
        self.path_stk = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/{}'.format(RunDir)
        self.sdb_stem = sdb_stem
        self.fpath = '{}/{}_SKS_SKKS.pairs'.format(self.path,self.sdb_stem)
        self.lam2 = [ ]
        #if kwargs are none:\

    def run(self):
        '''Function to build the pairs file'''
        # First match the SKS and SKKS
        start = ctime()
        print('Making Pairs')
        self.make_pairs()
        self.P = self.snr_check() # Foverwrite self.P eith only accepted events
        self.write_out(self.P)
        # Write initial pairs file so we can make piercepoints
        # Next generate the piercepoints and add them to the df
        print('Adding PiercePoints')
        self.gen_pp()
        self.add_pp()
        # Now calculate the d_SI vlaues
        print('Calculate difference in splitting intensity')
        self.add_DSI()
        # Finally stack the lamR surfaces to determine the lam2 values
        print('Calculate lambda 2 values')
        self.add_lam2()
        #Now test for matching and disrecpent pairs
        print('Apply 2-sigma test for discrepancy')
        self.match()
        # Apply a quick Signal to Noise test to get rid of the rreally bad data
        # print('{} pairs'.format(len(self.P)))

        # And save the result
        self.write_out(self.P)
        end = ctime()
        print('END. start {}, end {}'.format(start,end))

    def gen_pp(self):
        ''' Fucntion to test for whether the .pp file exists and if not call TauP to generate it and the corresponding mspp files '''
        if os.path.isfile('{}.pp'.format(self.fpath.split('.')[0])):
            #print(pf[:-6])
            self.pp = pd.read_csv('{}.pp'.format(self.fpath.split('.')[0]),delim_whitespace=True)
        else:
            print('Pierce Points file {}.pp doesnt not exist, calling pierce.sh'.format(self.fpath.split('.')[0]))
            p = call(shlex.split('/Users/ja17375/Shear_Wave_Splitting/Sheba/Programs/pierce.sh {}'.format(self.fpath)))
            self.pp = pd.read_csv('{}.pp'.format(self.fpath.split('.')[0]),delim_whitespace=True)
        # Load SDB and PP (pierce points data for a set of SKS-SKKS pairs)

        if os.path.isfile('{}.mspp'.format(self.fpath.split('.')[0])) is False:
            print('{}.mspp does not exist, creating'.format(self.fpath.strip('.pairs')))
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

    def write_out(self,df):
        df.to_csv('{}/{}_SKS_SKKS.pairs'.format(self.path,self.sdb_stem),sep=' ',index=False)

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
            # print('{}/{}/SKS/{}??_SKS.lam2'.format(self.path_stk,stat,fstem))
            if len(lam2_stem) is not 0:
                # I.e if glob has managed to find the sks lam2 surface file
                sks_lam2 = glob('{}/{}/SKS/{}??_SKS.{}'.format(self.path_stk,stat,fstem,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.{}'.format(self.path_stk,stat,fstem,ext))[0]
                Stk = Stacker(sks_lam2,skks_lam2,fstem,out)
                if mode == 'man':
                    self.lam2.append(Stk.lam2)
                elif mode == 'sheba':
                    self.lam2.append(Stk.sol[-1])
            else:
                fstem2 = '{}_{}'.format(stat,date)
                print('fstem2')
                # print('{}/{}/SKS/{}_*_SKS.{}'.format(self.path_stk,stat,fstem2,ext))
                sks_lam2 = glob('{}/{}/SKS/{}_*_SKS.{}'.format(self.path_stk,stat,fstem2,ext))[0]
                skks_lam2 = glob('{}/{}/SKKS/{}_*_SKKS.{}'.format(self.path_stk,stat,fstem2,ext))[0]
                # Now for a sanity check
                if (len(sks_lam2) is not 0) or (len(skks_lam2) is not 0):
                    Stk = Stacker(sks_lam2,skks_lam2,fstem2,out)
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

        match  = self.P[fast_test & lag_test] # Test for pairs that match within the given sigma range
        diff =  self.P.drop(index=match.index) # Remove matching pairs from original df to get the different pairs.
        # Write out matching and discepent dataframes
        match.to_csv('{}_matches.pairs'.format(fstem),index=False,sep=' ')
        diff.to_csv('{}_diffs.pairs'.format(fstem),index=False,sep=' ')
        # Open up mspp files
        mspp_match = open('{}_matches.mspp'.format(fstem),'w+')
        mspp_diff = open('{}_diffs.mspp'.format(fstem),'w+')

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

    def snr_check(self,t=2.0):
        ''' Run a quick and dirty check on Signal to Noise Ratio. If is is less than the threshold events are rejected'''

        self.accepted_i = [ ]
        self.d = [ ]
        print('There are {} pairs pre-SNR < 2 test'.format(len(self.P)))
        for i,row in self.P.iterrows():
            if row.SNR_SKS <= 2 or row.SNR_SKKS <= 2:
                #Test to see if Signal-to-Noise is too high
                print('SNR for SKS or SKKS less than 2, auto-reject')
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

    def __init__(self,df):
        self.df = df

    def plot_SNR(self):
        '''
        Make plots of SNR v dfast. In the style of phi_i v SNR from Restivo and Helffrich (2006)
        '''
        fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,sharey='row',figsize=(6,6))
        #Plot SNR for SKS
        ax1.plot(self.P.SNR_SKS,self.P.DFAST_SKS,'k.')
        ax1.set_ylabel('dfast')
        ax1.set_xlabel('S/N ratio')
        ax1.set_title('d fast SKS determination dependance on S/N')
        #Plot SNR for SKKS
        ax2.plot(self.P.SNR_SKKS,self.P.DFAST_SKKS,'k.')
        ax2.set_ylabel('dfast')
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
