#! /anaconda3/envs/splitwavepy/bin/python

################# stat_download.py ################
###################################################
# Script to download station data for temporary
# networks. Given a network code,
###################################################

# Just holding code i iused in ipython for nrows

stat,stla,stlo = [],[],[]

for i,stations in enumerate(s.networks[0].stations):
    ...:     sta = s.networks[0].stations[i]
    ...:     print('Stat: {}, STLA = {} STLO = {}'.format(sta.code,sta.latitude,sta.longitude))
    ...:     stat.append(sta.code)
    ...:     stla.append(sta.latitude)
    ...:     stlo.append(sta.longitude)

# Make pandas dataframe, double check station s are unique and pront 
