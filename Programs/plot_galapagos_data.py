#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 10:14:08 2021

@author: ja17375
"""
from obspy.taup import TauPyModel
import pygmt
import pandas as pd
from sphe_trig import vincenty_dist, vincenty_direct
from numpy import rad2deg
DATADIR='/Users/ja17375/DiscrePy/Data/GalapagosULVZ'
def calc_pp(evlo, evla, evdp, dist, azi):
    model = TauPyModel(model="iasp91")
    arr = model.get_pierce_points(evdp, dist, phase_list=['ScS'])
    for dep in arr[0].pierce:
        if dep[3] == 2889:
            pp_dist = rad2deg(dep[2])
            lat, lon = vincenty_direct(evla, evlo, azi, pp_dist)
            
    return lat, lon

if __name__ == '__main__':     
    region = [-115, -60, -30, 40]
    data = pd.read_csv(f'{DATADIR}/ScS/Galapagos_events_found.txt')
    model = TauPyModel(model="iasp91")
    
    fig = pygmt.Figure()
    fig.basemap(region=region, projection='M12c', frame='a10')
    fig.coast(shorelines='0.5p,black', resolution='l')
    for i, event in events.iterrows():
        for j, station in stations.iterrows():
            d, azi = vincenty_dist(data.EVLA, data.EVLO, station.stla, station.stlo)
            if d <= 60:
                ppla, pplo = calc_pp(data.EVLO, data.EVLA, data.EVDP, d, azi)
                # fig.plot(x=[data.EVLO, station.stlo], y=[data.EVLA, station.stla], pen='1p,black')
                fig.plot(x=pplo, y=ppla, style='c0.25c', color='blue', pen='0.5p,black')
                fig.plot(x=data.EVLO, y=.evla, style='a0.25c', color='yellow',pen='0.5p,black')
                fig.plot(x=station.stlo, y=station.stla, style='i0.2c', color='red', pen='0.5p,black')
    fig.show(method='external')
