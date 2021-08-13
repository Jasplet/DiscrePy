#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 11:50:05 2021

@author: ja17375
"""
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
from matplotlib.collections import PatchCollection
from matplotlib import colorbar
from matplotlib.patches import Polygon
import numpy as np
import matplotlib.pyplot as plt

DATA_DIR = '/Users/ja17375/DiscrePy/Sheba/Results/Combined/Filt_05Hz'

def global_heatmap(doms2plot, cmin=1):
    '''
    This is a function to make a heatmap showing our coverage. Each trigonal domain is coloured based off the number of phases that pass through it
    
    Args:
        geogeom_bins [array] - numpy array holding 
        
        cmin [int] - the minimum number of counts to draw a polygon for.
        
    Returns:
       fig [figure] - a 2 panel figure showing ScS and SnKS coverage in D'' and the Upper Mantle (reciever side)
    '''
    
    
    fig = plt.figure(figsize=(12,12))
    ax = fig.add_subplot(111,projection=ccrs.PlateCarree())
    # ax.set_extent(extent, crs=ccrs.PlateCarree())      
    ax.set_global()
    ax.add_feature(cfeature.GSHHSFeature(levels=[1],scale='high'))
    
    trig = []
    # fig = plt.gcf()
    for i, row in doms2plot.iterrows():
        vertices = make_vertices(row)
        trig += [Polygon(vertices,closed=True)]

    tps = PatchCollection(trig, alpha=0.6)
    tps.set_array(doms2plot.Count) # sets colors 
    tps.set_linestyle('-')
    tps.set_edgecolor('black')
    tps.set_linewidth(2.)
    ax.add_collection(tps) 
     # Add colorbar
    cax,kw =  colorbar.make_axes(ax,location='bottom',pad=0.05,shrink=0.7)
    cbar = fig.colorbar(tps, cax=cax, **kw)
    cbar.set_label('No. of phases in domain')
    ax.set_title(r"Coverage of SKS, SKKS phases in D$''$")

    plt.savefig('Global_phasecount_heatmap',format='png', dpi=500)

def make_vertices(dom):
    vertices = np.zeros([3,2])
    if (dom.V1_lon == 180) & (dom.lon < 0):
        V1_lon = -180
        if (dom.V2_lon == 180) & (dom.lon < 0):
            V2_lon = -180
            V3_lon = dom.V3_lon
        elif (dom.V3_lon == 180) & (dom.lon < 0):
            V2_lon = dom.V2_lon
            V3_lon = -180
        else:
            V2_lon = dom.V2_lon
            V3_lon = dom.V3_lon
    elif (dom.V2_lon == 180) & (dom.lon < 0):
        V1_lon = dom.V1_lon
        V2_lon = -180
        if (dom.V3_lon == 180) & (dom.lon < 0):
            V3_lon = -180
        else:
            V3_lon = dom.V3_lon
    elif (dom.V3_lon == 180) & (dom.lon < 0):
        V1_lon = dom.V1_lon
        V2_lon = dom.V2_lon
        V3_lon = -180
    else:
        V1_lon = dom.V1_lon
        V2_lon = dom.V2_lon
        V3_lon = dom.V3_lon
        
    vertices[:,0] = [V1_lon, V2_lon, V3_lon] #long
    vertices[:,1] = [dom.V1_lat, dom.V2_lat, dom.V3_lat]
    return vertices

if __name__ == '__main__':
    
    c_in = np.loadtxt('/Users/ja17375/DiscrePy/Sheba/Results/Combined/Combined_binned_all_T4.counts', skiprows=1)
    counts = c_in[:,0:2]
    cdf = pd.DataFrame(counts, columns=['bin_no','Count'])
    T4bins = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/Combined/binned_T4.tri', delim_whitespace=True)
    T4 = T4bins.drop_duplicates()
    domains2plot = pd.merge(cdf, T4, on='bin_no')
    global_heatmap(domains2plot)
    