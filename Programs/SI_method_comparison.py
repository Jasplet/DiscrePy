#! /usr/bin/env python
##############################
#   Program: SI_method_comparision.py
#
##############################
#   Author: J. Asplet
##############################
#   This program is designed to read sdb files that contain splitting analysis where
#   splitting intensity has been calculated by the low SNR approximation of Chevrot (2000) and
#   also by an implementation (made by J.Wookey) of the projection of radial derivatve onto the transverse
#   component (Chevrot 2000)
#
#   SI(Pr) = Splitting intensity by Projection
#   SI(Pa) = Splitting Inensity by Approximation

##############################
#   Import Statements
##############################
#   Standard Packages - all freely available
import os.path
import sys
import pandas as pd
import timeit
import shlex
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#######################################

sks = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SI_test/SI_by_proj_SKS_sheba_results.sdb',delim_whitespace=True)
skks = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SI_test/SI_by_proj_SKKS_sheba_results.sdb',delim_whitespace=True)

pairs = pd.read_csv('/Users/ja17375/Shear_Wave_Splitting/Sheba/Results/SI_test/SI_by_proj_05.pairs',delim_whitespace=True)

#def plot_si_Pr_v_Ap():
# '''Plot Splitting intensity calculated by approximatation against SI calculated by following the full projection method of Chevrot(20
# 00) '''
# Assume that sks  and skks file have been read in as pandas
fig1 = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(3,2)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[0,1])
ax3 = plt.subplot(gs[1,0])
ax4 = plt.subplot(gs[1,1])
ax5 = plt.subplot(gs[2,0])
ax6 = plt.subplot(gs[2,1])
####Plot SI for SPLIT data (Q > 0.7)
sks_split = sks[sks.Q > 0.7]
skks_split = skks[skks.Q > 0.7]
ax1.plot(sks_split['SI(Pa)'],sks_split['SI(Pr)'],'k.', label='sks')
ax1.plot(skks_split['SI(Pa)'],skks_split['SI(Pr)'],'r.',label='skks')
ax1.plot([-4, 4],[-4, 4],'k--',label='SI(Pr) = SI(Pa)')
ax1.set_ylim([-4,4])
ax1.set_xlim([-4,4])
ax1.set_title('Split Phases (Q > 0.7)')
ax1.set_xlabel(r'$\Delta SI$ by Approimxation')
ax1.set_ylabel(r'$\Delta SI$ by Projection')
ax1.legend()
## Beneath plot SI v Snr for both methods
ax3.plot(sks['Q'],sks['SI(Pr)'],'k.')

ax3.set_title('SI(Pr) v Q for SKS')
ax3.set_ylabel(r'$\Delta SI$ by Projection')
ax3.set_xlim([-1 ,1])
# Not for SI by proejction
ax5.plot(skks['Q'],skks['SI(Pr)'],'r.')
ax5.set_title('SI(Pr) v Q for SKKS')
ax5.set_ylabel(r'$\Delta SI$ by Projection')
ax5.set_xlim([-1 ,1])
#Plot SI for NULL data (Q < -0.7)
sks_null = sks[sks.Q < -0.7]
skks_null = skks[skks.Q < -0.7]
ax2.plot(sks_null['SI(Pa)'],sks_null['SI(Pr)'],'k.', label='sks')
ax2.plot(skks_null['SI(Pa)'],skks_null['SI(Pr)'],'r.',label='skks')
ax2.plot([-4, 4],[-2, 2],'k--',label='SI(Pr) = SI(Pa)')
ax2.set_ylim([-2,2])
ax2.set_xlim([-4,4])
ax2.set_xlabel(r'$\Delta SI$ by Approimxation')
ax2.set_ylabel(r'$\Delta SI$ by Projection')
ax2.set_title('Null Phases (Q < -0.7)')
ax2.legend()
## Beneath plot SI v Snr for both methods
ax4.plot(sks['Q'],sks['SI(Pa)'],'k.')
ax4.set_title('SI(Pa) v Q for SKS')
ax4.set_ylabel(r'$\Delta SI$ by Approimxation')
ax4.set_xlim([-1 ,1])
# Not for SI by proejction
ax6.plot(skks['Q'],skks['SI(Pa)'],'r.')
ax6.set_title('SI(Pa) v Q for SKKS')
ax6.set_ylabel(r'$\Delta SI$ by Approimxation')
ax6.set_xlim([-1 ,1])
fig2,ax = plt.subplots(1,1,figsize=(8,8))
# Make a plot of SI(Pr) v SI(Pa) for all data
ax.plot(sks['SI(Pa)'],sks['SI(Pr)'],'k.', label='sks')
ax.plot(skks['SI(Pa)'],skks['SI(Pr)'],'r.',label='skks')
ax.plot([-5, 5],[-5, 5],'k--',label='SI(Pr) = SI(Pa)')
ax.set_xlabel(r'$\Delta SI$ by Approimxation')
ax.set_ylabel(r'$\Delta SI$ by Projection')
ax.legend()
ax.set_ylim([-4,4])
ax.set_xlim([-4,4])
plt.show()

## Plot sks, skks Lambda 2 values and Delta SI values
fig3 = plt.figure(figsize=(10,10))
gs2 = gridspec.GridSpec(1,2)
ax1 = gs2[0]
ax2 = gs2[0]
# Plot Lam2 v dSI for sks
ax1.plot(pairs['LAM2_SKS'],pairs['SI(Pr)_SKS'],'k.',label = 'SI(Pr)')
ax1.plot(pairs['LAM2_SKS'],pairs['SI(Pa)_SKS'],'k.',label = 'SI(Pa)')
# Plot Lam2 v dSI for skks
ax2.plot(pairs['LAM2_SKKS'],pairs['SI(Pr)_SKKS'],'k.',label = 'SI(Pr)')
ax2.plot(pairs['LAM2_SKKS'],pairs['SI(Pa)_SKKS'],'k.',label = 'SI(Pa)')
#Add x/y labels
ax1.set_xlabel(r'$\lambda_{2}^{SKS}$')
ax2.set_xlabel(r'$\lambda_{2}^{SKKS}$')
ax1.set_ylabel(r'Splitting Intensity (SI)')
#Set ylimits
plt.show()
