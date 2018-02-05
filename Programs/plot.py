#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy as cart
import matplotlib.gridspec as gridspec
import obspy

def SKS_plot(file,title1):
    """
    Function to make diagnostice plots for a given file of splitting measuremtns
    """
    data = pd.read_csv(file,delim_whitespace=True)
    a = data['FAST']
    d = data.index[np.isnan(data['FAST']) == True].tolist() # Find any rows which contain NaNs
    data = data.drop(d)
    data = data[(data.QUAL != 'x')]
    fig,axs = plt.subplots(2, 2,sharex='col',figsize=(10,10))

    plt.subplot(221)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].FAST,yerr=data[(data.QUAL == 'n')].DFAST,fmt='kx',elinewidth=0.5,label='Null')
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].FAST,yerr=data[(data.QUAL != 'n')].DFAST,fmt='ko',elinewidth=0.5,label='Split')
    plt.legend(loc=2)

    plt.ylabel('Fast Direction (deg)')
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title('{} - Fast Direction'.format(title1))

    plt.subplot(223)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_FAST,yerr=data[(data.QUAL == 'n')].WL_DFAST,fmt='kx',elinewidth=0.5)
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_FAST,yerr=data[(data.QUAL != 'n')].WL_DFAST,fmt='ko',elinewidth=0.5)
    plt.ylim([-90,90])
    plt.yticks(np.arange(-90,91,30))
    plt.title('Jacks(Sheba) - Fast Direction')
    plt.xlabel('Back Azimuth')
    plt.ylabel('Fast Direction (deg)')

    plt.subplot(222)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].TLAG,yerr=data[(data.QUAL == 'n')].DTLAG,fmt='kx',elinewidth=0.5)
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].TLAG,yerr=data[(data.QUAL != 'n')].DTLAG,fmt='ko',elinewidth=0.5)
    plt.ylabel('Tlag (s)')
    plt.ylim([0,4])
    plt.title('{} - Lag Time'.format(title1))

    plt.subplot(224)
    plt.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_TLAG,yerr=data[(data.QUAL == 'n')].WL_DTLAG,fmt='kx',elinewidth=0.5)
    plt.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_TLAG,yerr=data[(data.QUAL != 'n')].WL_DTLAG,fmt='ko',elinewidth=0.5)
    plt.ylim([0,4])
    plt.ylabel('Tlag (s)')
    plt.xlabel('Back Azimuth')
    plt.title('Jacks(Sheba) - Lag Time')


    plt.tight_layout()
    plt.show()
