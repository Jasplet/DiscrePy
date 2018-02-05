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

    ax1 = subplot(221)
    ax1.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].FAST,yerr=data[(data.QUAL == 'n')].DFAST,fmt='kx',elinewidth=0.5,label='Null')
    ax1.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].FAST,yerr=data[(data.QUAL != 'n')].DFAST,fmt='ko',elinewidth=0.5,label='Split')
    ax1.legend(loc=2)

    ax1.ylabel('Fast Direction (deg)')
    ax1.ylim([-90,90])
    ax1.yticks(np.arange(-90,91,30))
    ax1.title('{} - Fast Direction'.format(title1))

    ax2 = subplot(223)
    ax2.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_FAST,yerr=data[(data.QUAL == 'n')].WL_DFAST,fmt='kx',elinewidth=0.5)
    ax2.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_FAST,yerr=data[(data.QUAL != 'n')].WL_DFAST,fmt='ko',elinewidth=0.5)
    ax2.ylim([-90,90])
    ax2.yticks(np.arange(-90,91,30))
    ax2.title('Jacks(Sheba) - Fast Direction')
    ax2.xlabel('Back Azimuth')
    ax2.ylabel('Fast Direction (deg)')

    ax3 = subplot(222)
    ax3.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].TLAG,yerr=data[(data.QUAL == 'n')].DTLAG,fmt='kx',elinewidth=0.5)
    ax3.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].TLAG,yerr=data[(data.QUAL != 'n')].DTLAG,fmt='ko',elinewidth=0.5)
    ax3.ylabel('Tlag (s)')
    ax3.ylim([0,4])
    ax3.title('{} - Lag Time'.format(title1))

    ax4 = subplot(224)
    ax4.errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_TLAG,yerr=data[(data.QUAL == 'n')].WL_DTLAG,fmt='kx',elinewidth=0.5)
    ax4.errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_TLAG,yerr=data[(data.QUAL != 'n')].WL_DTLAG,fmt='ko',elinewidth=0.5)
    ax4.ylim([0,4])
    ax4.ylabel('Tlag (s)')
    ax4.xlabel('Back Azimuth')
    ax4.title('Jacks(Sheba) - Lag Time')


    plt.tight_layout()
    plt.show()

def _lag(ax,baz,lag1,dlag1):

    ax.errorbar(baz,lag1,yerr=dlag1,fmt = 'kx',elinewidth=0.5)
    ax.ylim([0,4])
    ax.ylabel('Lag (s)')
    ax.xlabel('Back Azimuth (deg)')

def _fast(ax,baz,fast1,dfast1):

    ax.errorbar(baz,lag1,yerr=dlag1,fmt = 'kx',elinewidth=0.5)
    ax.ylim([0,4])
    ax.ylabel('Lag (s)')
    ax.xlabel('Back Azimuth (deg)')

def coverage(evla,evla,stla,stlo,stat):
    fig = plt.figure(figsize = (10,10))
    # We can either do and AzimuthalEquidistant projection centered on the staiton or a nice, Pacific-centered one
    # proj = cart.PlateCarree(central_longitude=180)
    proj = cart.AzimuthalEquidistant(central_longitude=stlo,central_latitude=stla)
    ax = plt.axes(projection=proj)
    ax.set_global() #This sets the axes extent to its maximum possible setting, so we can find these darn events

    ax.coastlines(resolution = '110m') # Coastline data is downloaded by Cartopy from Natural Earth (http://www.naturalearthdata.com)
    #Resolution options are 100m,50m and 10m

    ax.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')

    # Now add observed events
    ax.plot(evlo,evla,'ro',markersize = 5,transform = cart.Geodetic(),label='Event Location')

    # ax.set_xticks([-130,-125,-120,-115,-110,-105,-100], crs=proj)
    # ax.set_yticks([30,35,40,45,50,55,60], crs=proj)
    stat = 'NEW'
    ax.plot(stlo,stla,'kv',transform=cart.Geodetic(),markersize=10,label='Station Loction')
    plt.title('Coverage for Station {}'.format(stat))
    ax.legend()
