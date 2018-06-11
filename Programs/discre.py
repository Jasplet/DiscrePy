#! /anaconda3/envs/splitwavepy/bin/python

# Welcome to discre.py. This script (or maybe module??) is for testing whether
# SKS SKKS pairs exhibit matching or discrepent splitting. This will be done
# by using the splitting measured for each phase to correct the splitting on
# the other.
### Imports
import numpy as np
import pandas as pd
import obspy as ob
#import splitwavepy as sw
import matplotlib.pyplot as plt
from stack import Stacker,plot_stack
from glob import glob
from datetime import datetime
# Maybe some others

###############################################################################

def pair_stack(pairs,path):
    ''' Runs Stacker for all the desired pairs (a .pairs file)'''
    lam2 = []
    print('Running')
    for i,f in enumerate(pairs.DATE.values):
        #rint('It {}, time is {} '.format(i,str(datetime.now())))
        # First get the right DATE,TIME and STATION
        date,time,stat = pairs.DATE[i], pairs.TIME[i], pairs.STAT[i]
        fstem = '{}_{}_{}'.format(stat,date,time)

        lam2_stem = glob('{}/{}/SKS/{}??_SKS.lam2'.format(path,stat,fstem))
        print('{}/{}/SKS/{}??_SKS.lam2'.format(path,stat,fstem))
        if len(lam2_stem) is not 0:
            # I.e if glob has managed to find the sks lam2 surface file
            sks_lam2 = glob('{}/{}/SKS/{}??_SKS.lam2'.format(path,stat,fstem))[0]
            skks_lam2 = glob('{}/{}/SKKS/{}??_SKKS.lam2'.format(path,stat,fstem))[0]
            Stk = Stacker(sks_lam2,skks_lam2,fstem)
            lam2.append(Stk.sol[-1])
        else:
            fstem2 = '{}_{}'.format(stat,date)
            sks_lam2 = glob('{}/{}_*_SKS.lam2'.format(path,fstem2))[0]
            skks_lam2 = glob('{}/{}_*_SKS.lam2'.format(path,fstem2))[0]
            # Now for a sanity check
            if (len(sks_lam2) is not 0) or (len(skks_lam2) is not 0):
                Stk = Stacker(sks_lam2,skks_lam2)
                lam2.append(Stk.sol[-1])
            else:
                #print('lam2 surfaces cannot be found, skipping')
                pass
#        Now lets get the lambda 2 values
    print('Lam2 max: {} Lam2 min: {}'.format(max(lam2),min(lam2)))
    return lam2

def plot_lam2(x,lam2):
    print('Plotting')

    plt.plot(x,lam2,'k.')
    plt.ylabel('lambda 2 values')
    plt.yticks(np.arange(0,2,step=0.2))
    plt.ylim([0,2])
    plt.show()

def write_lam2(pairs,lam2):
    '''Adds lam2 values to pairs'''
    l2df = {'LAM2' : lam2}
    ldf = pd.DataFrame(l2df)
    pairs['LAM2'] = ldf

    pairs.to_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Accepted_SKS_SKKS_all_w_lam2.pairs',sep=' ')
    return pairs

def show_stacks(df):
    '''Function to find and plot desired surface stacks based on the LAMDA2 value '''
    ### Plot Min Lamnda 2
    paths = [ ]
    p_stem = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split'
    minl2 = df[df.LAM2 == df.LAM2.min()]
    date,time,stat = minl2.DATE.values[0], minl2.TIME.values[0], minl2.STAT.values[0]
    l_stem = '{}_{}_{}'.format(stat,date,time)
    paths.append('{}/{}/{}'.format(p_stem,stat,l_stem))

    maxl2 = df[df.LAM2 == df.LAM2.max()]
    date,time,stat = maxl2.DATE.values[0], maxl2.TIME.values[0], maxl2.STAT.values[0]
    l_stem = '{}_{}_{}'.format(stat,date,time)
    paths.append('{}/{}/{}'.format(p_stem,stat,l_stem))

    plot_stack(paths)


if __name__ == '__main__':
    print('This is discre.py')
    date_time_convert = {'TIME': lambda x: str(x),'DATE': lambda x : str(x)}
    p = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/Jacks_Split/Accepted_SKS_SKKS_all.pairs',delim_whitespace=True,converters=date_time_convert)
    path = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Runs/Jacks_Split'
    lam2 =  pair_stack(p,path)
    p2 = write_lam2(p,lam2) # p2 contians lam2 values

    show_stacks(p2)
    #print(lam2)
    #plot_lam2(p.index.values,lam2)
