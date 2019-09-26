#! /usr/bin/env python

#path_plot.py
####################
# Script to plot S/ScS and SKS/SKKS raypath through the earth
##################################################################
#Read in S, ScS
import numpy as np
import matplotlib.pyplot as plt

theta_S,r_S = np.loadtxt('S.path',skiprows=1,unpack=True)
theta_ScS,r_ScS = np.loadtxt('ScS.path',skiprows=1,unpack=True)
# Read in SKS, SKKS
theta_SKS,r_SKS = np.loadtxt('SKS.path',skiprows=1,unpack=True)
theta_SKKS,r_SKKS = np.loadtxt('SKKS.path',skiprows=1,unpack=True)


ax = plt.subplot(111,projection='polar')
ax.set_theta_zero_location('N')
# Plot some basic layers we know exist in the earth
theta = np.radians(np.linspace(0,360,200))
moho = 6350*np.ones(200)
l_440 = (6371-440)*np.ones(200)
l_660 = (6371-660)*np.ones(200)
cmb = 3481*np.ones(200)
d_prime = 3731*np.ones(200)
inner_core = 1217*np.ones(200)
#ax.plot(theta,moho,c='darkgrey',linestyle='dotted')
ax.plot(theta,l_440,c='darkgrey',linestyle='dotted')
ax.plot(theta,l_660,c='darkgrey',linestyle='dotted')
ax.plot(theta,d_prime,c='darkgrey',linestyle='dashed')
ax.plot(theta,cmb,c='grey',linewidth=2)
ax.plot(theta,inner_core,c='darkgrey',linestyle='dotted')
# Now plot SKS and SKKS paths (or S / ScS)
ax.plot(np.radians(theta_S),r_S,c='blue')
ax.plot(np.radians(theta_ScS),r_ScS,c='blue')
# Add SKS/sKKS
ax.plot(np.radians(theta_SKS-105),r_SKS,c='red')
ax.plot(np.radians(theta_SKKS-105),r_SKKS,c='red')
ax.set_rmax(6371)
ax.grid(False)
ax.set_yticklabels([])
ax.tick_params(labelcolor='white')
# ax.set_thetamax(90)
# ax.set_thetamin(30)
plt.savefig('SKS_SKKS_S_ScS_wedgepaths.png',transparent=True,format='png',dpi=400)
plt.show()
