#! /usr/bin/env python

"""
This is a script to test the basic concepts of building interative functionality for a plot window
"""

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)

#for n in range(1,5):
ax.plot([i**2 for i in range(0,10)])

def nextplot(event):
    if event.key == 'a':
        print ("This is a good plot")
        return event.key
        plt.close()
    elif event.key == 'x':
        print ("NEXT")
        fig.canvas.mpl_disconnect(cid)
        plt.close()

cid = fig.canvas.mpl_connect('key_press_event',nextplot)
print(cid)
plt.show()


#fig.canvas.mpl_disconnect(cid)
#
# import numpy as np
# import matplotlib.pyplot as plt
#
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title('click on points')
#
# line, = ax.plot(np.random.rand(100), 'o', picker=5)  # 5 points tolerance
#
# def onpick(event):
#     thisline = event.artist
#     xdata = thisline.get_xdata()
#     ydata = thisline.get_ydata()
#     ind = event.ind
#     points = tuple(zip(xdata[ind], ydata[ind]))
#     print('onpick points:', points)
#
# fig.canvas.mpl_connect('pick_event', onpick)
#
# plt.show()
