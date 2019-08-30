import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os
import sys


def animate(i):
    ax1.clear()
    for root, dirs, files in os.walk("./plotter"):
        files = [ fi for fi in files if fi.endswith(".txt") ]
        for filename in files:
            filename = "./plotter/" + filename
            pullData = open(filename,"r").read()
            dataArray = pullData.split('\n')
            xar = []
            yar = []
            for eachLine in dataArray:
                if len(eachLine)>1:
                    x,y = eachLine.split(',')
                    xar.append(float(x))
                    yar.append(float(y))
            ax1.plot(xar,yar)

if __name__ == '__main__':
	output = sys.argv[1]
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)

	animate(1)
	fig.savefig("plots/"+output)

