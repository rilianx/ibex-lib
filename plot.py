import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import math
import os
import random
from stat import *
import multiprocessing
from tkinter import *
import time
from itertools import islice

#this file was created to plot particles from PSO algorithm implemented into IBEX
# python3 plot.py

def main():

    #Create the base plot
    plot()

    #Call a function that reads a file and draws on plot
    plotdata()

    plt.show()
    print('Done')

def plot():
	global line,ax1,canvas,fig

	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)
	line, = ax1.plot([], [])

#format must be 2 dimension vector, variables separated by (,) and every vector separated by (;)
#i.e. 1,2;4,5;78.324324,-12342 every line is a iterations from PSO
def plotdata():
	with open('output.txt') as f:
		for line in f:
			Px = []
			Py = []
			bestx = []
			besty = []
			aux = line.replace('inf', 'math.inf') #change math format
			aux = aux.replace('nan', '0') #change math format
			aux = aux.replace('\n', '') #delete symbol
			particulas = aux.split(';') #separate particles positions
			print(particulas)
			for particula in particulas[:-1]:
				# for each particle's position, separate x and y						
				s = particula.split(",")
				Px.append(float(s[0]))
				Py.append(float(s[1]))
			sb = particulas[-1].split(",")
			bestx.append(float(sb[0]))
			besty.append(float(sb[1]))

			fig.clear() #clear plot for next iteration
			ax1.plot()
			plt.plot(Px,Py, 'r.', markersize=1)
			plt.plot(bestx,besty, 'b.', markersize=2) #best in blue
			plt.pause(0.5) #wait	
	f.close()
if __name__ == '__main__':
    main()
