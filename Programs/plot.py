#! /usr/bin/env python
### Script containing varous plotting functions for splitting Measurements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as cart
import cartopy
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


    axs[0,0].errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].FAST,yerr=data[(data.QUAL == 'n')].DFAST,fmt='kx',elinewidth=0.5,label='Null')
    axs[0,0].errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].FAST,yerr=data[(data.QUAL != 'n')].DFAST,fmt='ko',elinewidth=0.5,label='Split')
    axs[0,0].legend(loc=2)

    axs[0,0].set_ylabel('Fast Direction (deg)')
    axs[0,0].set_ylim([-90,90])
    axs[0,0].set_yticks(np.arange(-90,91,30))
    axs[0,0].set_title('{} - Fast Direction'.format(title1))

    axs[0,1].errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_FAST,yerr=data[(data.QUAL == 'n')].WL_DFAST,fmt='kx',elinewidth=0.5)
    axs[0,1].errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_FAST,yerr=data[(data.QUAL != 'n')].WL_DFAST,fmt='ko',elinewidth=0.5)
    axs[0,1].set_ylim([-90,90])
    axs[0,1].set_yticks(np.arange(-90,91,30))
    axs[0,1].set_title('Jacks(Sheba) - Fast Direction')
    axs[0,1].set_xlabel('Back Azimuth')
    axs[0,1].set_ylabel('Fast Direction (deg)')

    _lag(axs[1,0],data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].TLAG,data[(data.QUAL == 'n')].DTLAG,fmt='kx')
    _lag(axs[1,0],data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].TLAG,data[(data.QUAL != 'n')].DTLAG,fmt='ko')
    # axs[1,0].set_ylabel('Tlag (s)')
    # axs[1,0].set_ylim([0,4])
    # axs[1,0].set_title('{} - Lag Time'.format(title1))

    axs[1,1].errorbar(data[(data.QUAL == 'n')].BAZ,data[(data.QUAL == 'n')].WL_TLAG,yerr=data[(data.QUAL == 'n')].WL_DTLAG,fmt='kx',elinewidth=0.5)
    axs[1,1].errorbar(data[(data.QUAL != 'n')].BAZ,data[(data.QUAL != 'n')].WL_TLAG,yerr=data[(data.QUAL != 'n')].WL_DTLAG,fmt='ko',elinewidth=0.5)
    axs[1,1].set_ylim([0,4])
    axs[1,1].set_ylabel('Tlag (s)')
    axs[1,1].set_xlabel('Back Azimuth')
    axs[1,1].set_title('Jacks(Sheba) - Lag Time')

    plt.tight_layout()
    plt.show()



def SKKS_plot(stat,phase):

    data = pd.read_csv('/Users/ja17375/Python/Shear_Wave_Splitting/Measurements/{}_{}_Splitting.txt'.format(stat,phase),delim_whitespace=True)
    a = data['FAST']
    d = data.index[np.isnan(data['FAST']) == True].tolist() # Find any rows which contain NaNs
    data = data.drop(d)
    data = data[(data.QUAL != 'x')]
    fig = plt.figure(figsize = [20,20])
    axs = []
    gs = gridspec.GridSpec(2,3)
    proj = cart.AzimuthalEquidistant(central_longitude=data.STLO[0],central_latitude=data.STLA[0])


    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[:,1:],projection = proj)

    ax1.errorbar(data.BAZ,data.TLAG,yerr = data.DTLAG,fmt='kx',elinewidth=0.5)
    ax1.set_ylim([0,4])
    ax1.set_ylabel('Lag (s)')
    ax1.set_xlabel('Back Azimuth (deg)')
    # _lag(ax1,data.BAZ,data.TLAG,data.DTLAG,'kx')
    # _fast(ax2,data.BAZ,data.FAST,data.DFAST,'kx')
    # coverage(ax3,data.EVLA,data.EVLO,data.STLA[0],data.STLO[0],stat)
    # plt.show()
    ax2.errorbar(data.BAZ,data.FAST,yerr=data.DFAST,fmt='kx',elinewidth=0.5)
    ax2.set_ylim([-90,90])
    ax2.set_ylabel('Fast Direction (s)')
    ax2.set_xlabel('Back Azimuth (deg)')
    ### Plot Coverage on third axis object

    ax3.set_global() #This sets the axes extent to its maximum possible setting, so we can find these darn events
    ax3.coastlines(resolution = '110m') # Coastline data is downloaded by Cartopy from Natural Earth (http://www.naturalearthdata.com)
    #Resolution options are 100m,50m and 10m
    ax3.add_feature(cartopy.feature.OCEAN, zorder=0)
    ax3.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    # Now add observed events
    ax3.plot(data.EVLO,data.EVLA,'ro',markersize = 5,transform = cart.Geodetic(),label='Event Location')
    # ax.set_xticks([-130,-125,-120,-115,-110,-105,-100], crs=proj)
    # ax.set_yticks([30,35,40,45,50,55,60], crs=proj)
    ax3.plot(data.STLO,data.STLA,'kv',transform=cart.Geodetic(),markersize=10,label='Station Loction')
    ax3.set_title('Coverage for Station {}'.format(stat))
    ax3.legend()
    plt.show()

def plot_lag(ax,baz,lag1,dlag1,fmt):

    ax.errorbar(baz,lag1,yerr=dlag1,fmt=fmt,elinewidth=0.5)

    ax.set_ylim([0,4])
    ax.set_ylabel('Lag (s)')
    ax.set_xlabel('Back Azimuth (deg)')

def plot_fast(ax,baz,fast1,dfast1,fmt):

    ax.errorbar(baz,fast1,yerr=dfast1,fmt = fmt,elinewidth=0.5)
    ax.set_ylim([-90,90])
    ax.set_ylabel('Fast Direction (s)')
    ax.set_xlabel('Back Azimuth (deg)')

def coverage(ax,evla,evlo,stla,stlo,stat):

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
    ax.set_title('Coverage for Station {}'.format(stat))
    ax.legend()
