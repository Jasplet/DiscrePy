#! /usr/bin/env python
##############################
#   Program: sheba_plot.py
#
##############################
#   Author: J. Asplet
##############################
#   A python program to plot Shear_Wave_Splitting measurements
#   made using Sheba. This program assumes that the results files
#   have already been cleaned up into final_resuts and STAT_results
#   textfiles using tidyup_by_stat.sh and tidyup.sh
##############################
#   Import Statements
##############################
#   Standard Packages - all freely available
import obspy as ob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as cart
import cartopy
import matplotlib.gridspec as gridspec

import time
##############################

class plotter:

    def __init__(self,filename):
        #Load data into object on initiation
        self.path = '/Users/ja17375/Shear_Wave_Splitting/Sheba/Results'
        self.data = pd.read_csv('{}/{}'.format(self.path,filename),delim_whitespace=True)


    def plot_by_stat(self,station,save=False):
        """
        Function to plot by station IF final results file has been read in
        """
        fig = plt.figure()
        axs = []
        gs = gridspec.GridSpec(2,4) # Creates a 2x3 grid in the figure space for plotting in
        data = self.data[self.data.STAT=='station']
        proj = cart.AzimuthalEquidistant(central_longitude=self.data.STLO[0],central_latitude=self.data.STLA[0])

        ax1 = plt.subplot(gs[0,0:2])
        ax2 = plt.subplot(gs[1,0:2])
        ax3 = fig.add_subplot(gs[:,2:],projection = proj) # have to use add_subplots in order to add a different projection

        ax1.errorbar(self.data.BAZ,self.data.FAST,yerr=self.data.DFAST,fmt = 'bx',elinewidth=0.5,label='Fast Direction')
        ax1.set_ylim([-90,90])
        ax1.set_ylabel('Fast Direction (deg)')
        ax1.set_xlabel('Back Azimuth (deg)')
        ax1.legend()

        ax2.errorbar(self.data.BAZ,self.data.TLAG,yerr=self.data.DTLAG,fmt = 'bx',elinewidth=0.5,label='Lag ')
        ax2.set_ylim([0,4])
        ax2.set_ylabel('Lag (s)')
        ax2.set_xlabel('Back Azimuth (deg)')
        ### Plot Coverage on third axis object
        #ax3.set_extent([-180,180,-60,90])
        ax3.set_global() #This sets the axes extent to its maximum possible setting, so we can find these darn events
        ax3.coastlines(resolution = '110m') # Coastline data is downloaded by Cartopy from Natural Earth (http://www.naturalearthdata.com)
        #Resolution options are 100m,50m and 10m
        #ax3.add_feature(cartopy.feature.OCEAN, zorder=0)
        #ax3.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
        # Now add observed events

        ax3.plot(self.data.EVLO,self.data.EVLO,'bx',markersize = 5,transform = cart.Geodetic(),label='Event Locations')
        # ax.set_xticks([-130,-125,-120,-115,-110,-105,-100], crs=proj)
        # ax.set_yticks([30,35,40,45,50,55,60], crs=proj)
        ax3.plot(self.data.STLO[0],self.data.STLA[0],'kv',transform=cart.Geodetic(),markersize=10,label='Station {}'.format(stat))
        ax3.set_title('SKS coverage for Station {}'.format(stat))
        ax3.legend()


        if save is True:
            pass
            # plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Python/Figures/{}_SKS_SKKS_plot'.format(stat))
            # print('Saving')
        else:
            plt.show()
