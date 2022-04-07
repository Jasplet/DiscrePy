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

sks_raw = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/SI_proj/SI_proj_SKS_sheba_results.sdb',delim_whitespace=True)
skks_raw = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/SI_proj/SI_proj_SKKS_sheba_results.sdb',delim_whitespace=True)

sks = sks_raw[sks_raw.SNR > 5]
skks = skks_raw[skks_raw.SNR > 5]

#def plot_si_Pr_v_Ap():
# '''Plot Splitting intensity calculated by approximatation against SI calculated by following the full projection method of Chevrot(20
# 00) '''
# Assume that sks  and skks file have been read in as pandas
fig1 = plt.figure(figsize=(10,14))
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
ax1.set_xlabel(r'$SI$ by Approximation')
ax1.set_ylabel(r'$SI$ by Projection')
ax1.legend()
## Beneath plot SI v Snr for both methods
ax3.plot(sks['Q'],sks['SI(Pr)'],'k.')

ax3.set_title('SKS')
ax3.set_ylabel(r'$SI$ by Projection')
ax3.set_xlim([-1 ,1])
# Not for SI by proejction
ax5.plot(skks['Q'],skks['SI(Pr)'],'r.')
ax5.set_title('SKKS')
ax5.set_ylabel(r'$SI$ by Projection')
ax5.set_xlim([-1 ,1])
ax5.set_xlabel(r'Q')
ax6.set_xlabel(r'Q')
#Plot SI for NULL data (Q < -0.7)
sks_null = sks[sks.Q < -0.7]
skks_null = skks[skks.Q < -0.7]
ax2.plot(sks_null['SI(Pa)'],sks_null['SI(Pr)'],'k.', label='sks')
ax2.plot(skks_null['SI(Pa)'],skks_null['SI(Pr)'],'r.',label='skks')
ax2.plot([-4, 4],[-4, 4],'k--',label='SI(Pr) = SI(Pa)')
ax2.set_ylim([-2,2])
ax2.set_xlim([-4,4])
ax2.set_xlabel(r'$SI$ by Approximation')
ax2.set_ylabel(r'$SI$ by Projection')
ax2.set_title('Null Phases (Q < -0.7)')
ax2.legend()
## Beneath plot SI v Snr for both methods
ax4.plot(sks['Q'],sks['SI(Pa)'],'k.')
ax4.set_title('SKS')
ax4.set_ylabel(r'$SI$ by Approximation')
ax4.set_xlim([-1 ,1])
# Not for SI by proejction
ax6.plot(skks['Q'],skks['SI(Pa)'],'r.')
ax6.set_title('SKKS')
ax6.set_ylabel(r'$SI$ by Approximation')
ax6.set_xlim([-1 ,1])

# Adjust space between figures so there is enough room
fig1.subplots_adjust(hspace=0.25)
# Save figure properly
fig1.savefig('/Users/ja17375/DiscrePy/Figures/SI_Pr_v_Pa_6panel.eps',dpi = 400, format = 'png')
######## END FIG1    ###########

######## START FIG2  ###########
fig2,ax = plt.subplots(1,1,figsize=(8,8))
# Make a plot of SI(Pr) v SI(Pa) for all data
ax.plot(sks['SI(Pa)'],sks['SI(Pr)'],'k.', label='sks')
ax.plot(skks['SI(Pa)'],skks['SI(Pr)'],'r.',label='skks')
ax.plot([-5, 5],[-5, 5],'k--',label='SI(Pr) = SI(Pa)')
ax.set_xlabel(r'$\Delta SI$ by Approximation')
ax.set_ylabel(r'$\Delta SI$ by Projection')
ax.legend()
ax.set_ylim([-4,4])
ax.set_xlim([-4,4])
# Save figure
fig2.savefig('/Users/ja17375/DiscrePy/Figures/SI_Pr_v_Pa.png',dpi = 400, format = 'png')
######## END FIG 2   ##########

######## Start Fig 3 ##########
pairs = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/E_pacific/E_pacific_05.pairs',delim_whitespace=True)
diff = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/E_pacific/E_pacific_05_diffs_l2.pairs',delim_whitespace=True)
match = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/E_pacific/E_pacific_05_matches_l2.pairs',delim_whitespace=True)
nulls = pd.read_csv('/Users/ja17375/DiscrePy/Sheba/Results/E_pacific/E_pacific_05_nulls.pairs',delim_whitespace=True)
## Plot sks, skks Lambda 2 values and Delta SI values
fig3 = plt.figure(figsize=(8,16))
gs = gridspec.GridSpec(3,2)
#Make bottom row first so we can share x axes without gettingn undefined errors
ax3 = fig3.add_subplot(gs[1,0])
ax4 = fig3.add_subplot(gs[1,1],sharey=ax3)

ax1 = fig3.add_subplot(gs[0,0],sharex=ax3)
ax2 = fig3.add_subplot(gs[0,1],sharey=ax1,sharex=ax4)

ax5 = fig3.add_subplot(gs[2,0])
ax6 = fig3.add_subplot(gs[2,1],sharey=ax5)
# Plot Lam2 v SI, dSI for sks (top row)
#SI by approximation
# ax1.plot(pairs['LAM2_SKS'],pairs['SI(Pa)_SKS'],'k.',label = 'SI(Pa) SKS')
ax1.plot(nulls['LAM2_SKS'],nulls['SI(Pa)_SKS'],'k.',label = 'discrepant')
# ax1.plot(diff['LAM2_SKS'],diff['SI(Pa)_SKS'],'b.',label = 'discrepant')
# ax1.plot(match['LAM2_SKS'],match['SI(Pa)_SKS'],'r.',label = 'matching' )
# ax1.legend()

#SI by projection
# ax2.plot(pairs['LAM2_SKS'],pairs['SI(Pr)_SKS'],'k.',label = 'SI(Pr) SKS')
ax2.plot(nulls['LAM2_SKS'],nulls['SI(Pr)_SKS'],'k.',label = 'SI(Pr) SKS')
# ax2.plot(diff['LAM2_SKS'],diff['SI(Pr)_SKS'],'b.',label = 'SI(Pr) SKS')
# ax2.plot(match['LAM2_SKS'],match['SI(Pr)_SKS'],'r.',label = 'SI(Pr) SKS')
# Plot Lam2 v SI,dSI fro SKKS (bottom row)
#SI by approximation
# ax3.plot(pairs['LAM2_SKKS'],pairs['SI(Pa)_SKKS'],'k.',label = 'SI(Pa) SKKS')
ax3.plot(nulls['LAM2_SKKS'],nulls['SI(Pa)_SKKS'],'k.',label = 'SI(Pa) SKKS')
# ax3.plot(diff['LAM2_SKKS'],diff['SI(Pa)_SKKS'],'b.',label = 'SI(Pa) SKKS')
# ax3.plot(match['LAM2_SKKS'],match['SI(Pa)_SKKS'],'r.',label = 'SI(Pa) SKKS')
#SI by projection
# ax4.plot(pairs['LAM2_SKKS'],pairs['SI(Pr)_SKKS'],'k.',label = 'SI(Pr) SKKS')
ax4.plot(nulls['LAM2_SKKS'],nulls['SI(Pr)_SKKS'],'k.',label = 'SI(Pr) SKKS')
# ax4.plot(diff['LAM2_SKKS'],diff['SI(Pr)_SKKS'],'b.',label = 'SI(Pr) SKKS')
# ax4.plot(match['LAM2_SKKS'],match['SI(Pr)_SKKS'],'r.',label = 'SI(Pr) SKKS')
# L2bar v dSI (approx)
# ax5.plot(pairs['LAM2_BAR'],pairs['D_SI_Pa'],'k.', label = r'$\Delta SI (Pa)')
ax5.plot(nulls['LAM2_BAR'],nulls['D_SI_Pa'],'k.', label = r'$\Delta SI (Pa)')
# ax5.plot(diff['LAM2_BAR'],diff['D_SI_Pa'],'b.', label = r'$\Delta SI (Pa)')
# ax5.plot(match['LAM2_BAR'],match['D_SI_Pa'],'r.', label = r'$\Delta SI (Pa)')
ax5.plot([0,1],[0.4,0.4],'k--')
# L2bar v dSI (proj)
# ax6.plot(pairs['LAM2_BAR'],pairs['D_SI_Pr'],'k.', label = r'$\Delta SI (Pr)')
ax6.plot(nulls['LAM2_BAR'],nulls['D_SI_Pr'],'k.', label = r'$\Delta SI (Pr)')
# ax6.plot(diff['LAM2_BAR'],diff['D_SI_Pr'],'b.', label = r'$\Delta SI (Pr)')
# ax6.plot(match['LAM2_BAR'],match['D_SI_Pr'],'r.', label = r'$\Delta SI (Pr)')
ax6.plot([0,1],[0.4,0.4],'k--')
#Add x/y labels
ax1.set_ylabel(r'Splitting Intensity (SI)')

ax3.set_xlabel(r'$\lambda_{2}$')
ax3.set_ylabel(r'Splitting Intensity (SI)')
ax4.set_xlabel(r'$\lambda_{2}$')
ax5.set_ylabel(r'$\Delta SI$')
ax5.set_xlabel(r'$\bar{\lambda_{2}}$')
ax6.set_xlabel(r'$\bar{\lambda_{2}}$')
# ax6.set_ylabel(r'$\Delta SI$')
#Add titles
ax1.set_title(r'SI by approximation for SKS')
ax2.set_title(r'SI by projection for SKS')
ax3.set_title('SI by approximation for SKKS')
ax4.set_title('SI by projection for SKKS')
ax5.set_title(r'$\Delta SI (Pa) V \bar{\lambda_2}$')
ax6.set_title(r'$\Delta SI (Pr) V \bar{\lambda_2}$')
#Set ylimits and xlimits
ax1.set_ylim([-2,2.5])
# ax2.set_ylim([-2,2.5])
ax3.set_ylim([-3,3])
ax3.set_xlim([0,0.14]) # 0.14 for all pairs
ax4.set_xlim([0,0.14]) # 0.08 for splits
ax4.set_ylim([-3,3])
ax5.set_ylim([0,3]) #3.5 for all pairs, 1 for splits
ax6.set_ylim([0,3]) #2.5 for all pairs, 1 for splits
ax5.set_xlim([0,0.25]) #0.25 for all pairs
ax6.set_xlim([0,0.25])

fig3.subplots_adjust(hspace=0.25)
fig3.savefig('/Users/ja17375/DiscrePy/Figures/SI_comp_null_pairs.png',dpi=400,format='png')
########## END FIG3 ###########

######### START FIG4 ##########
fig4,(ax1,ax2) = plt.subplots(2,1,figsize=(8,8))
# Make a plot of SI(Pr) v SI(Pa) for all data
Pa_Pr_diff_sks = abs(sks['SI(Pa)']  - sks['SI(Pr)'])
Pa_Pr_diff_skks = abs(skks['SI(Pa)'] - skks['SI(Pr)'])
bins = np.arange(0,6,0.25)
ax1.set_title(r'Histograms of |SI(Approximation) - SI(Projection)| for SKS and SKKS')
ax1.hist(Pa_Pr_diff_sks,bins)
ax2.hist(Pa_Pr_diff_skks,bins)

ax1.set_xlabel(r'$|SI(Pa) - SI(Pr)|$ SKS')
ax2.set_xlabel(r'$|SI(Pa) - SI(Pr)|$ SKKS')
# ax.set_ylim([,4])
# ax.set_xlim([-4,4])
fig4.savefig('/Users/ja17375/DiscrePy/Figures/SI_Pr_Pa_difference_histograms.png',dpi=400,format='png')
plt.show()
