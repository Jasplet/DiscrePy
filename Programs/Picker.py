import numpy as np
import obspy
import matplotlib.pyplot as plt

class WindowPicker:
    """
    Pick a Window
    """

    def __init__(self,st,wbeg1,wbeg2,wend1,wend2,t0,**kwargs):

        #t0 = 60 seconds before traveltime (aka the start of the trimming seismogram)
        self.st = st # Obspy stream containing BHN and BHE
        self.tt = t0 + 60  # traveltime (predicted)
        self.delta = st[0].stats.delta
        self.t = [t0 + self.delta*i for i,v in enumerate(self.st[0].data)]
        # make initial window rnages attributes
        self.wbeg1 = wbeg1
        self.wbeg2 = wbeg2
        self.wend1 = wend1
        self.wend2 = wend2
        # Base plot (before interactive stuff)
        fig,self.ax = plt.subplots(figsize = (12,8))
        self.canvas = fig.canvas
        # Add seismograms
        self.ax.plot(self.t,self.st[0].data,label='BHN')
        self.ax.plot(self.t,self.st[1].data,label='BHE')
        # window limit lines
        self.wbeg1line = self.ax.axvline(self.wbeg1,linewidth=1,color='r',visible=True)
        self.wbeg2line = self.ax.axvline(self.wbeg2,linewidth=1,color='r',visible=True)
        self.wend1line = self.ax.axvline(self.wend1,linewidth=1,color='r',visible=True)
        self.wend2line = self.ax.axvline(self.wend2,linewidth=1,color='r',visible=True)
        self.cursorline= self.ax.axvline(100,linewidth=1,color='0.5',visible=False)
        self.pred_tt= self.ax.axvline(self.tt,linewidth=1,color='k',visible=True)
        _,self.ydat = self.wbeg1line.get_data()
        # set limits
        lim_max = max([self.st[0].data.max(),self.st[1].data.max()]) * 1.1
        lim_min = min([self.st[0].data.min(),self.st[1].data.min()])* 1.1
        # self.ax.set_aspect('equal')
        self.ax.set_ylim([lim_min,lim_max])
        self.ax.set_xlim(t0,max(self.t)) # Set ylim in relative time (from stsrt of stream )

        self.connect() # Dev only
        plt.show()

    def connect(self):
        # self.cidclick = self.canvas.mpl_connect('button_press_event', self.click)
        self.cidmotion = self.canvas.mpl_connect('motion_notify_event', self.motion)
        # self.cidrelease = self.canvas.mpl_connect('button_release_event', self.release)
        self.cidenter = self.canvas.mpl_connect('axes_enter_event', self.enter)
        self.cidleave = self.canvas.mpl_connect('axes_leave_event', self.leave)
        self.cidkey = self.canvas.mpl_connect('key_press_event', self.keypress)

    # def click(self,event):
    #     if event.inaxes is not self.ax: return
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
        '''
        if event.key == "a":
            print('WBEG Start')
            self.x1 = event.xdata
            self.wbeg1line.set_data(self.x1,self.ydat)
            self.canvas.draw()
            print(self.x1)
        if event.key == "d":
            print('WBEG End')
            self.x2 = event.xdata
            self.wbeg2line.set_data(self.x2,self.ydat)
            self.canvas.draw()
            print(self.x2)
        if event.key == "z":
            print('WEND Start')
            self.x3 = event.xdata
            self.wend1line.set_data(self.x3,self.ydat)
            self.canvas.draw()
            print(self.x3)
        if event.key == "c":
            print('WEND End')
            self.x4 = event.xdata
            self.wend2line.set_data(self.x4,self.ydat)
            self.canvas.draw()
            print(self.x4)
        if event.key == "q":
            print('DISCONNECT')
            self.disconnect()

    def enter(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.cursorline.set_visible(True)
        self.canvas.draw()

    def leave(self,event):
        if event.inaxes is not self.ax: return
        self.cursorline.set_visible(False)
        self.canvas.draw()

    def motion(self,event):
        if event.inaxes is not self.ax: return
        x = event.xdata
        self.cursorline.set_data([x,x],self.ydat)
        self.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        # self.canvas.mpl_disconnect(self.cidclick)
        self.canvas.mpl_disconnect(self.cidmotion)
        self.canvas.mpl_disconnect(self.cidenter)
        self.canvas.mpl_disconnect(self.cidleave)
        self.wbeg1, self.wbeg2, self.wend1,self.wend2 = sorted((self.x1, self.x2,self.x3,self.x4))
        plt.close()
