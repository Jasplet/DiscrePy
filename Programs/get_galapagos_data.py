#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 16:41:26 2021
Downloads sac data from IRIS using event IDs
Assumes event and station metadata has been downloaded/found usign Wilbur (or similar)

@author: ja17375
"""
import os
import pandas as pd
import obspy 
from obspy.clients.fdsn.header import (FDSNNoDataException)
from obspy.taup import TauPyModel
from sphe_trig import vincenty_dist
DATA_DIR='/Users/ja17375/DiscrePy/Data/GalapagosULVZ'


def download_waveforms(event, station):
    
    client = obspy.clients.fdsn.client.Client('IRIS') 
    dist, azi = vincenty_dist(event.evla, event.evlo, station.stla, station.stlo)
    if (dist >= 60) & (dist <= 85):
        # Start by only looking for ScS between 60 and 85 deg
        scs_tt = predict_scs_tt(event.evdp, dist)
        evttime = obspy.core.utcdatetime.UTCDateTime(event.Time)
        start = evttime + scs_tt - 180
        end = evttime + scs_tt + 180
        try:
            st = client.get_waveforms(network=station.Net,
                                      station=station.Stat,
                                      channel='BHE,BHN,BHZ',
                                      location='*',
                                      starttime=start,
                                      endtime=end,
                                      attach_response=True)
            if len(st) == 3:
                for trace in st:
                    save_trace(trace, evttime, station.Stat)
                return True
            else:
                print('Unexpected number of traces')
                print(st)
                print(len(st))
        except FDSNNoDataException:
            print(f'No Data for {station.Stat}')
        
    else:
        print(f'Event-Station distance {dist:4.3f} is too great or too small.')
        return False
def predict_tt(evdp, dist, phase='ScS'):
    model = TauPyModel()
    # uses IASP91 by default
    arrivals = model.get_travel_times(source_depth_in_km=evdp,
                                      distance_in_degree=dist,
                                      phase_list=[phase])
    tt = arrivals[0].time
    return tt

def save_trace(tr, time, stat):

    try:
        os.mkdir(f'{DATA_DIR}/ScS/{stat}')
    except FileExistsError:
        print(f'{DATA_DIR}/{stat} already exists, Hooray! Less work for me!')
    tr_id = f'{stat}_{time.year:04d}{time.julday:03d}_{time.hour:02d}{time.minute:02d}{time.second:02d}.{tr.stats.channel}'
    tr.write(f'{DATA_DIR}/ScS/{stat}/{tr_id}',format='SAC', byteorder=1)


if __name__ == '__main__':  
    events = pd.read_csv(f'{DATA_DIR}/Deep_SouthAmerica_events.txt')
    stations = pd.read_csv(f'{DATA_DIR}/ScS_stations.txt')
    outfile = open(f'{DATA_DIR}/ScS/Galapagos_events_found.txt','w')
    outfile.write('EventID,Time,Mag,MagType,EVLA,EVLO,EVDP,Station,Network,STLA,STLO\n')
    for i, event in events.iterrows():
        if (event.magnitude >= 5.7) & (event.magnitude <= 6.7):
            print(f'{event.Time} EVLA:{event.evla} EVLO:{event.evlo}, Mag:{event.magnitude}')
            for j, station in stations.iterrows():
                data_found = download_waveforms(event, station)
                if data_found:
                    outfile.write('{},{},{},{},{:5.4f},{:5.4f},{:5.4f},{},{},{:5.4f},{:5.4f}\n'.format(event.EventID, event.Time,
                                                                         event.magnitude, event.MagType, event.evla,
                                                                         event.evlo, event.evdp, station.Stat, station.Net,
                                                                         station.stla, station.stlo
                                                                         ))
    outfile.close()
    assert outfile.closed