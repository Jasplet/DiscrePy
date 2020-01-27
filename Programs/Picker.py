import numpy as np
import obspy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.fftpack import fft

class WindowPicker:
    """
    Pick a Window start/end range, for use with cluster analysis code
    """

    def __init__(self,st,wbeg1,wbeg2,wend1,wend2,tt,**kwargs):

        #t0 = 60 seconds before traveltime (aka the start of the trimming seismogram)
        self.st = st # Obspy stream containing BHN and BHE
        # st.plot()
        self.tt = tt
        t0 = tt - 60
        self.delta = st[0].stats.delta
        self.t = [t0 + self.delta*i for i,v in enumerate(self.st[0].data)]
        # make initial window rnages attributes
        self.wbeg1 = wbeg1
        self.wbeg2 = wbeg2
        self.wend1 = wend1
        self.wend2 = wend2
        # Base plot (before interactive stuff)
        fig = plt.figure(figsize = (18,8))
        gs = gridspec.GridSpec(2,2)
        self.ax1 = plt.subplot(gs[1,:]) # Bottown Row, for window picking
        self.ax2 = plt.subplot(gs[0,:]) # Top Row, for fft plot
        self.plot_fft()
        self.canvas = fig.canvas
        # Add seismograms
        self.ax1.plot(self.t,self.st[0].data,label='BHN',color='darkorange')
        self.ax1.plot(self.t,self.st[1].data,label='BHE',color='dodgerblue')
        self.ax1.set_xlabel('Time relative to origin (s)')
        # Add legend
        self.ax1.legend()
        # window limit lines
        self.wbeg1line = self.ax1.axvline(self.wbeg1,linewidth=2,color='r',visible=True)
        self.wbeg2line = self.ax1.axvline(self.wbeg2,linewidth=2,color='r',visible=True)
        self.wend1line = self.ax1.axvline(self.wend1,linewidth=2,color='g',visible=True)
        self.wend2line = self.ax1.axvline(self.wend2,linewidth=2,color='g',visible=True)
        self.cursorline= self.ax1.axvline(100,linewidth=1,color='0.5',visible=False)
        self.pred_tt= self.ax1.axvline(self.tt,linewidth=1,color='k',visible=True)
        _,self.ydat = self.wbeg1line.get_data()

        # set limits
        self.lim_max = max([self.st[0].data.max(),self.st[1].data.max()]) * 1.1
        self.lim_min = min([self.st[0].data.min(),self.st[1].data.min()])* 1.1
        # self.ax1.set_aspect('equal')
        self.ax1.set_ylim([self.lim_min,self.lim_max])
        self.ax1.set_xlim(t0,max(self.t)) # Set ylim in relative time (from stsrt of stream )
        # Add some labels
        self.phaselabel = self.ax1.text(self.tt + 1,self.lim_max*0.8,"IASP91\nPred.\nArrival",multialignment='left')
        self.wbeg1label = self.ax1.text(self.wbeg1 - 3, self.lim_min*0.85,'S',color='r',fontsize=14)
        self.wbeg2label = self.ax1.text(self.wbeg2 - 3, self.lim_min*0.85,'F',color='r',fontsize=14)
        self.wend1label = self.ax1.text(self.wend1 - 3, self.lim_min*0.85,'S',color='g',fontsize=14)
        self.wend2label = self.ax1.text(self.wend2 - 3, self.lim_min*0.85,'F',color='g',fontsize=14)
        print("'a' & 'd' set the window beginnning range")
        print("'z' & 'c' set the window end range")
        self.connect() # Dev only
        plt.show()

    def plot_fft(self):
        ''' Takes and fft of both components (within wbeg1 and wend2 and plots). This does not get updated by changes in windowing'''
        # Trim traces to intial wbeg and wend. Use delta to work out the indexes that correspond to the window start/end positions.
        # The transform these positions to integers so we can use slicing
        st = self.st.copy()

        bhe = st[0]
        bhn = st[1]

        # bhe.trim(bhe[0].stats.starttime + 45,bhe[0].stats.starttime + 90)
        # bhn.trim(bhn[0].stats.starttime + 45,bhn[0].stats.starttime + 90)

        #bhe  = self.st[1].data[n_start:n_end] # Slice out the data in the window
        N = len(bhn) # number of samples in traces
        # Set sample spacing in f-domain
        df = 1.0/(2.0*self.delta) # Sample frequency
        xf = np.linspace(0.0,df,N/2) # Frequencies up to F_nyquist (N/2*df)
        # Take the fft
        BHN_f = fft(bhn.data)
        BHE_f = fft(bhe.data)
        # print(N//2)
        # print(BHN_f[1:N//2].shape)
        # print(xf[1:N//2].shape)
        #Now plot the spectra
        self.ax2.semilogy(xf[1:N//2], 2.0/N * np.abs(BHN_f[1:N//2]), color='darkorange')
        self.ax2.semilogy(xf[1:N//2], 2.0/N * np.abs(BHE_f[1:N//2]), color='dodgerblue')
        self.ax2.legend(['BHN', 'BHE'])
        self.ax2.set_xlabel('Frequency [Hz]')
        self.ax2.set_xlim([0,2.0])
        plt.grid()
        # plt.show()

    def connect(self):
        # self.cidclick = self.canvas.mpl_connect('button_press_event', self.click)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion)
        # self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        self.cidenter = self.canvas.mpl_connect('axes_enter_event', self.enter)
        self.cidleave = self.canvas.mpl_connect('axes_leave_event', self.leave)
        self.cidkey = self.canvas.mpl_connect('key_press_event', self.keypress)

    # def click(self,event):
    #     if event.inaxes is not self.ax1: return
    #     x = event.xdata
    #     if event.button == 1:
    #         self.wbeg1 = x
    #         self.wbeg1line.set_data([x,x],self.ydat)
    #         self.canvas.draw()
    #     if event.button == 3:
    #         self.wbeg2 = x
    #         self.wbeg2line.set_data([x,x], self.ydat)
    #         self.canvas.draw()

    def keypress(self,event):
        ''' Define a set of keypress responses
        'a' & 'd' set the window beginnning range
        'z' & 'c' set the window end range
        'q' exit the plot and returns the current WBEG, WEND
        The vertical line markers and annotations are redrawn after each Key Press
        '''
        # print("'a' & 'd' set the window beginnning range")
        # print("'z' & 'c' set the window end range")
        if event.key == "a":
            print('WBEG Start')
            self.x1 = event.xdata
            self.wbeg1line.set_data(self.x1,self.ydat)
            self.wbeg1label.set_position((self.x1 - 3, self.lim_min*0.85))
            self.canvas.draw()
            print(self.x1)
        if event.key == "d":
            print('WBEG End')
            self.x2 = event.xdata
            self.wbeg2line.set_data(self.x2,self.ydat)
            self.wbeg2label.set_position((self.x2 - 3, self.lim_min*0.85))
            self.canvas.draw()
            print(self.x2)
        if event.key == "z":
            print('WEND Start')
            self.x3 = event.xdata
            self.wend1line.set_data(self.x3,self.ydat)
            self.wend1label.set_position((self.x3 - 3, self.lim_min*0.85))
            self.canvas.draw()
            print(self.x3)
        if event.key == "c":
            print('WEND End')
            self.x4 = event.xdata
            self.wend2line.set_data(self.x4,self.ydat)
            self.wend2label.set_position((self.x4 - 3, self.lim_min*0.85))
            self.canvas.draw()
            print(self.x4)
        if event.key == "w":
            print('Bad/noisey waveform, will tell SHEBA to skip')
            self.x1 = False
            self.disconnect()
        if event.key == "q":
            print('DISCONNECT')
            self.disconnect()

    def enter(self,event):
        if event.inaxes is not self.ax1: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.cursorline.set_visible(True)
        self.canvas.draw()

    def leave(self,event):
        if event.inaxes is not self.ax1: return
        self.cursorline.set_visible(False)
        self.canvas.draw()

    def motion(self,event):
        if event.inaxes is not self.ax1: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        # self.canvas.mpl_disconnect(self.cidclick)
        self.canvas.mpl_disconnect(self.cidmotion)
        self.canvas.mpl_disconnect(self.cidenter)
        self.canvas.mpl_disconnect(self.cidleave)
        if self.x1 == False:
            'Bad waveform that we want to skip'
            self.wbeg1, self.wbeg2, self.wend1,self.wend2 = (None,None,None,None)
        else:
            self.wbeg1, self.wbeg2, self.wend1,self.wend2 = sorted((self.x1, self.x2,self.x3,self.x4))
        plt.close()
