import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation
import os
import math
import random
import subprocess
import multiprocessing
import traceback
import argparse


# No estoy usando esto
# =========================
class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class Box:
    def __init__(self, box_string):
        box_string = box_string.replace('inf', "math.inf")
        box_string = box_string.replace('nan', "0")
        box_var = eval(box_string)
        self.createBox(box_var['pts'][0], box_var['pts'][1],
                       box_var['diam_x'], box_var['diam_y'])

    def __repr__(self):
        return str((self.x, self.y))

    def createBox(self, x, y,  diam_x, diam_y):
        self.x = x
        self.y = y
        self.diam_x = diam_x
        self.diam_y = diam_y


class DataDict:
    box_set = []
    particles_set = []
    particle_pbest = []
    particle_gbest = []
    ub_set = []
# =========================


def plot():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # ax.set_aspect(1)
    line, = plt.plot([], [], 'b-', label='particle')
    line2, = plt.plot([], [], 'r-', label='pbest')
    plt.title('Plotter')
    plt.xlim(-300, 300)
    plt.ylim(-300, 300)
    return fig, line, ax, line2


def updateplot(num, q, l, ax, l2):
    try:
        result = q.get_nowait()
        if result != 'Q':
            particle_x = []
            particle_y = []
            pbest_x = []
            pbest_y = []
            for p in [patches.Rectangle(
                        (box['pts'][0], box['pts'][1]),
                        box['diam_x'], box['diam_y'],
                        fill=False,
                        edgecolor='black',
                        linestyle='solid',
                        lw=0.1
                        ) for box in result[0]]:
                ax.add_patch(p)
            for p in result[1]:
                particle_x.append(p[0])
                particle_y.append(p[1])
            for p in result[2]:
                pbest_x.append(p[0])
                pbest_y.append(p[1])
            gbest_x = [result[3][0]]
            gbest_y = [result[3][1]]
            ub_x = [result[4][0]]
            ub_y = [result[4][1]]
            l.set_data(particle_x, particle_y)
            line_gbest = plt.Line2D(gbest_x, gbest_y, marker='.'
                                    lw=0.5, markeredgecolor='black')
            line_ub = plt.Line2D(ub_x, ub_y, marker='.',
                                 lw=0.2, markeredgecolor="green")
            l2.set_data(pbest_x, pbest_y)
            ax.add_line(line_gbest)
            ax.add_line(line_ub)
            return tuple(ax) + (l, l2)
        else:
            q.close()
    except Exception as e:
        # print("Esto es un error:")
        # print(e)
        return tuple(patch_list) + (l, ) + (l2, ) + tuple(line_list)


def readFile():
    data = DataDict()
    f = open("output.txt")
    reader = f.read()
    reader = reader.replace('inf', "math.inf")
    reader = reader.replace('nan', "0")
    reader = reader.split("\n")
    for it in reader:
        q.put((
            eval(it[0]),  # box
            eval(it[1]),  # particle
            eval(it[2]),  # pbest
            eval(it[3]),  # gbest
            eval(it[4]),  # ub
        ))


def main():
    try:
        q = multiprocessing.Queue()
        fig, line, ax, line2 = plot()
        # Read the file
        readFile()
        line_ani = animation.FuncAnimation(
            fig, updateplot, q.qsize(),
            fargs=(q, line, ax, line2),
            interval=100, blit=True, repeat=False
        )
        plt.show()
        print("Done")
    except Exception as e:
        print(e)
        print("Hay un error")


if __name__ == '__main__':
    main()
