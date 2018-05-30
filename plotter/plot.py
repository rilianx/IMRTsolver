import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

def animate(i):
    ax1.clear()
    for root, dirs, files in os.walk("."):
        files = [ fi for fi in files if fi.endswith(".txt") ]
        for filename in files:
            pullData = open(filename,"r").read()
            dataArray = pullData.split('\n')
            xar = []
            yar = []
            for eachLine in dataArray:
                if len(eachLine)>1:
                    x,y = eachLine.split(',')
                    xar.append(int(x))
                    yar.append(int(y))
            ax1.plot(xar,yar)

ani = animation.FuncAnimation(fig, animate, interval=1000)
plt.show()
