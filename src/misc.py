#! /usr/bin/env python
######################
# misc.py
######################
# A library of misc functions that are useful, but may be required in several programs. Therefore I have dencentralised the here.
# Contents:
#   std_norm - a function to calculate the value of the STANDARD normal distribution for a given x value
#   gen_norm - function to evaluate a generalised normal distribution for a given x
#   trunc_norm - evalutes a trucated normal distribution

#####################################################################################################
# Imports
import math
import scipy
import numpy as np
import matplotlib.pyplot as plt

def std_norm(x):
    ''' Evalutate the standard normal distibution (mean =0 , var = 1) for a given x'''
    reg = (1/math.sqrt(2*math.pi))
    e = math.exp(((-1.)*(x ** 2)/2.))
    z = reg*e

    return z

def gen_norm(u,s,x):
     ''' Evaluate the general normal distibution (mean =u , std = s) for a given x'''
     reg = (1/(s*math.sqrt(2*math.pi)))
     e = math.exp(((-1)*(x - u)**2)/2*(s**2))
     z = reg*e
     return z

def cdf(xi):
    ''' Evalute a given Cumulative Denssity function for x '''
    c = 0.5*(1 + math.erf(xi/math.sqrt(2)))
    return c

def trunc_norm(u,s,x,a=0,b=np.inf):
    ''' Evaluate a truncated normal distribution (mean = u, var = v, lower bound = a, upper bound = b) for a given x'''

    if a < x < b:
        # If x is within truncation limits
        xi = (x- u)/s
        xi_a = (a -u)/s
        xi_b = (b -u)/s
        num = std_norm(xi) # Call general normal function
        dem = cdf(xi_b) - cdf(xi_a) # Evalutate differences in CFDs at truncation bounds.
        # print(num, dem)
        z = num/dem
    else:
        z = 0
    return z

def plot_trunc(u,s,x,a=0,b=np.inf):
    ''' Plots a given truncated normal distribution'''
    d_trunc = [ ] # List to hold the distirbution
    d_gen = [ ]
    for i in x: #Iterate over a range of X values
        d_trunc.append(trunc_norm(u,s,a,b,x=i))
        d_gen.append(gen_norm(u,s,x=i))
    plt.plot(x,d_gen,'--')
    plt.plot(x,d_trunc,color='darkorange',linestyle='dashed')
    plt.xlabel(r'$X$')
    plt.ylabel(r'$\Psi(x)$')
    plt.title('Truncated normal')
    plt.show()

def plot_windows(st,st2):
    ''' FUnction to plot out the windowing logic hat we use for a given pair
    st = obpsy stream object trimmed for SKS & filtered
    sr2 = obspy sream object trimmed for SKKS & filtered
    '''
 # Assume that we already have the SKS and SKKS cut traces read in as st and st2
    #Only want N and E components.
    # We already have the traveltimes calculated. From how I cut the traces I know that the traveltime pick is 60 seconds from the start
    # SO we plot the SKS trace and then need to add on SKKS
    # Diff in the traveltimes is: (user1 hold traveltime accoring to TauP using IASPEI91)
    shift = st2[0].stats.sac.user1 - st[0].stats.sac.user1
    print('TT diff between SKS and SKKS is',shift)
    rel_skks = 60 + shift
    fig = plt.figure(figsize = (12,8))
    st.plot(fig=fig,type='relative')
    (ax1,ax2) = fig.axes[0],fig.axes[1]
    # For SKS
    ax1.plot([45,45],[-3000,3000],'k--')
    ax1.plot([60,60],[-3000,3000],'k-')
    ax1.plot([75,75],[-3000,3000],'k--')
    ax1.plot([90,90],[-3000,3000],'k--')
    #ax1.set_ylim([-2000,3000])
    #For SKKS
    ax1.plot([rel_skks,rel_skks],[-3000,3000],'r-')
    ax1.plot([rel_skks-15,rel_skks-15],[-3000,3000],'r--')
    ax1.plot([rel_skks+15,rel_skks+15],[-3000,3000],'r--')
    ax1.plot([rel_skks+30,rel_skks+30],[-3000,3000],'r--')
    ax1.set_ylim([st[1].data.min(),st[1].data.max()])
    ax1.tick_params(axis='both',which='major',labelsize= 15)

    ax2.plot([45,45],[-4000,3000],'k--')
    ax2.plot([60,60],[-4000,3000],'k-')
    ax2.plot([75,75],[-4000,3000],'k--')
    ax2.plot([90,90],[-4000,3000],'k--')
    #ax1.set_ylim([-2000,3000])
    ax2.plot([rel_skks,rel_skks],[-4000,3000],'r-')
    ax2.plot([rel_skks-15,rel_skks-15],[-4000,3000],'r--')
    ax2.plot([rel_skks+15,rel_skks+15],[-4000,3000],'r--')
    ax2.plot([rel_skks+30,rel_skks+30],[-4000,3000],'r--')
    ax2.set_ylim([st[0].data.min(),st[0].data.max()])
    ax2.tick_params(axis='both',which='major',labelsize= 15)
    plt.savefig('/Users/ja17375/Shear_Wave_Splitting/Figures/Raw_Figs/Window_Logic_V1.eps',format='eps',dpi=400)
    plt.show()
