#!/usr/bin/env python
import sys
import splitwavepy as sw
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

def plotall(filestem):

    plt.figure(figsize=(15,10))
    gs = gridspec.GridSpec(7, 6)
    ax0 = plt.subplot(gs[0, 0:5])
    plt.ion()
    plt.show()

    # plot trace at top (this will put sks windows in automatically)
    sks = sw.load(filestem + '_sks.trnm')
    sks.srcpoldata()._ptr(ax0)
    # add the skks window
    skks = sw.load(filestem + '_skks.trnm')
    ax0.axvline(skks.data.wbeg(), linewidth=1, color='k')
    ax0.axvline(skks.data.wend(), linewidth=1, color='k')

     # plot FFT power spectrum from radial component
    ax06 = plt.subplot(gs[0,5])
    power0 = np.abs(np.fft.fft(sks.srcpoldata().data()[0]))**2
    power1 = np.abs(np.fft.fft(sks.srcpoldata().data()[1]))**2
    freq = np.fft.fftfreq(sks.srcpoldata().data().shape[1], sks.srcpoldata().delta)
    ax06.fill_between(freq, power0, label='rad', alpha=0.5)
    ax06.fill_between(freq, power1, label='trans', alpha=0.5)
    ax06.set_xlim(0, 0.6)
    ax06.set_xlabel('Frequency (Hz)')
    ax06.set_ylabel('Power')
    ax06.legend()

    # now start plotting particle motions and error surf
    ax10 = plt.subplot(gs[1,0])
    sks.data._ppm(ax10)
    ax20 = plt.subplot(gs[2,0])
    sks.data_corr()._ppm(ax20)
    ax11 = plt.subplot(gs[1:3, 1:3])
    sks._psurf(ax11, vals=sks.energy1/sks.energy2, conf95=True, marker=True, info=True)
    # ax10.axis('off')
    # ax20.axis('off')
    # ax11.axis('off')

    ax13 = plt.subplot(gs[1,3])
    skks.data._ppm(ax13)
    ax23 = plt.subplot(gs[2,3])
    skks.data_corr()._ppm(ax23)
    ax14 = plt.subplot(gs[1:3, 4:6])
    skks._psurf(ax14, vals=skks.energy1/skks.energy2, conf95=True, marker=True, info=True)

    # now load eigenm files
    sks = sw.load(filestem + '_sks.eigm')
    skks = sw.load(filestem + '_skks.eigm')
    # now start plotting particle motions and error surf
    ax30 = plt.subplot(gs[3,0])
    sks.data._ppm(ax30)
    ax40 = plt.subplot(gs[4,0])
    sks.data_corr()._ppm(ax40)
    ax31 = plt.subplot(gs[3:5, 1:3])
    sks._psurf(ax31, vals=sks.lam1/sks.lam2, conf95=True, marker=True, info=True)

    ax33 = plt.subplot(gs[3,3])
    skks.data._ppm(ax33)
    ax43 = plt.subplot(gs[4,3])
    skks.data_corr()._ppm(ax43)
    ax34 = plt.subplot(gs[3:5, 4:6])
    skks._psurf(ax34, vals=skks.lam1/skks.lam2, conf95=True, marker=True, info=True)

    # now load crossm files
    sks = sw.load(filestem + '_sks.xcrm')
    skks = sw.load(filestem + '_skks.xcrm')
    # now start plotting particle motions and error surf
    ax50 = plt.subplot(gs[5,0])
    sks.data._ppm(ax50)
    ax60 = plt.subplot(gs[6,0])
    sks.data_corr()._ppm(ax60)
    ax51 = plt.subplot(gs[5:7, 1:3])
    sks._psurf(ax51, vals=sks.fisher(), conf95=True, marker=True, info=True)

    ax53 = plt.subplot(gs[5,3])
    skks.data._ppm(ax53)
    ax63 = plt.subplot(gs[6,3])
    skks.data_corr()._ppm(ax63)
    ax54 = plt.subplot(gs[5:7, 4:6])
    skks._psurf(ax54, vals=skks.fisher(), conf95=True, marker=True, info=True)


    # plt.tick_params(
    #     axis='both',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom='off',      # ticks along the bottom edge are off
    #     top='off',         # ticks along the top edge are off
    #     labelbottom='off')

    plt.tight_layout()
    plt.draw()
    qual = input('Data Quality ([g]ood/[o]k/[p]oor) \n >')
    plt.draw()
    disc  = input('Are SKS/SKKS discrepant (can add comments)? ([y]es/[n]o/[u]ndeterminable) \n ?')
    plt.close()


    return qual,disc

if __name__ == '__main__':
    # print(sys.argv[1])

    plotall(sys.argv[1])
    plt.show()
