from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from itertools import islice
from tkinter import *
from stat import *

import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import multiprocessing
import argparse
import random
import time
import math
import os

#this file was created to plot particles from PSO algorithm implemented into IBEX
#python3 plot.py -f <file_name.txt> -x_1 <int> -x_2 <int>
#text must be formated like x_1,x_2,...,x_n;x_1,x_2,...,x_n and end line for each iteration.
#file must be indicated on execution.
#will only plot 2 var at a time, you must indicate each one on execution.

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
	parser = argparse.ArgumentParser()
	parser.add_argument("--file", "-f", type=str, required=True)
	parser.add_argument("--x_1", "-x_1", type=str, required=True)
	parser.add_argument("--x_2", "-x_2", type=str, required=True)
	args = parser.parse_args()

	print("file: "+args.file)
	print("x_1: "+args.x_1)
	print("x_2: "+args.x_2)
	x_1 = int(args.x_1)
	x_2 = int(args.x_2)
	with open(args.file) as f:
		for line in f:
			try:
				input("Press any key")
			except SyntaxError:
				pass
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
				# for each particle's position, separate x_1 and x_2
				s = particula.split(",")
				Px.append(float(s[x_1]))
				Py.append(float(s[x_2]))
				sb = particulas[-1].split(",")
				bestx.append(float(sb[x_1]))
				besty.append(float(sb[x_2]))

			fig.clear() #clear plot for next iteration
			ax1.plot()
			plt.ylim(-5,5)
			plt.xlim(-5,5)
			#plt.plot([-5,5],[-3,7], 'y-', markersize=0.5) # x - y < -2
			#plt.plot([-5,5],[-2,8]) # x - y > -3
			plt.plot(Px,Py, 'r.', markersize=1)
			plt.plot(bestx,besty, 'b.', markersize=2) #best in blue
	f.close()
if __name__ == '__main__':
    main()
