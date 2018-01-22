#!/usr/bin/env python

#from obspy.core import read
from obspy.clients.iris import Client
from obspy import core

stat = 'AAM' # Station Code
net = 'US' # Station Network
ch = ['BHZ','BH1','BH2'] # Channel
loc = '00' # Location 
# These 4 parameter identify the trace to download


start = core.UTCDateTime("2017-09-19T18:14:38")
end = core.UTCDateTime( start + 1800)

client = Client() # Not sure what this does but its gets rid of the "self" error message
str = core.stream.Stream()
for x in range(0,3):
 str += client.timeseries(network=net,station=stat,channel=ch[x],location=loc,starttime=start,endtime=end) 
## Plot Traces
start_plot = start + 222 #Manipulate starttime to get a nicer plot
end_plot = start_plot + 960
str.plot(size = [800,600],starttime = start_plot,endtime = end_plot, type = 'relative')