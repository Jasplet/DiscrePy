#! /anaconda3/envs/splitwavepy/bin/python

# Welcome to discre.py. This script (or maybe module??) is for testing whether
# SKS SKKS pairs exhibit matching or discrepent splitting. This will be done
# by using the splitting measured for each phase to correct the splitting on
# the other.
### Imports
import numpy as np
import pandas as pd
import obspy as ob
import splitwavepy as sw
import matplotlib.pyplot as plt
from stack import Stacker
from glob import glob
# Maybe some others

###############################################################################
if __name__ == ' __main__':

    p = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Resuts/Jacks_Split/Accepted_SKS_SKKS_all.pairs')
    path = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split/'
    pair_stack(p,path)
    # call function


def pair_stack(pairs,path):
    ''' Runs Stacker for all the desired pairs (a .pairs file)'''

    for i,f in enumerate(pairs.DATE.values):
        print('It {}, time is {} '.format(i,str(datetime.now())))
        # First get the right DATE,TIME and STATION
        date,time,stat = pairs.DATE[i], pairs.TIME[i], pairs.STAT[i]
        fstem = '{}_{}_{}'.format(stat,date,time)

        sks_lam2 = glob('{}/{}/SKS/{}??_SKS.lam2'.format(path,stat,fstem))
        if len(sks_lam2) is not 0:
            # I.e if glob has managed to find the sks lam2 surface file
            skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.lam2'.format(path,stat,fstem))
            Stk = Stacker(sks_lam2,skks_lam2)
        else:
            fstem2 = '{}_{}'.format(stat,date)
            sks_lam2 = glob('{}/{}_*_SKS.lam2'.format(path,fstem2))
            skks_lam2 = glob('{}/{}_*_SKS.lam2'.format(path,fstem2))
            # Now for a sanity check
            if (len(sks_lam2) is not 0) or (len(skks_lam2) is not 0):

            else:
                print('lam2 surfaces cannot be found, skipping')
                pass
