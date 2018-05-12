import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import random as rand
import time as t
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import scipy.integrate as scint

steps = 1000
finalt = 10
N = 6
m = 1
b = 0.1 #damping
g = 10 # gravity
k = 40 # drag force constant
time = np.linspace(0, finalt, steps+1)
dt = finalt/steps

""" Class definition """

# the cloth class is defined here
# this doesn't create the variables so the definitions can't be defined on other variables in here
class clothPoint:
    def __init__(self,row,col):
        self.m = 1 # mass
        self.left = False # single out edge cases
        self.top = False
        self.bottom = False
        self.right = False

        # position from inputs
        self.row = row
        self.col = col

        self.xpos = np.zeros(steps)
        self.ypos = np.zeros(steps)

# pointlist contains objects of type clothPoint
pointlist = []
for i in range(N):
    for j in range(N):
        pointlist.append(clothPoint(i,j))
pointlist = np.reshape(pointlist,(N,N)) # reshapes pointlist to NxN

# define edge cases and fill in the booleans from clothPoint for each object
for i in range(N):
    for j in range(N):
        if i == 0:
            pointlist[i][j].top = True
        if i == (N-1):
            pointlist[i][j].bottom = True
        if j == 0:
            pointlist[i][j].left = True
        if j == (N-1):
            pointlist[i][j].right = True

# initial conditions
for i in range(N):
    for j in range(N):
        pointlist[i][j].xpos[0] = -i # initial x position
        pointlist[i][j].ypos[0] = -i # initial y position

""" Runga Kutta method to solve for new positions """

def rk4(h,xn,yn,f,i,j,k):
    k1 = h*f(xn,yn,i,j,k)
    k2 = h*f(xn+h/2,yn+k1/2,i,j,k)
    k3 = h*f(xn+h/2,yn+k2/2,i,j,k)
    k4 = h*f(xn+h,yn+k3,i,j,k)

    return yn+k1/6+k2/3+k3/3+k4/6

# first derivative of x
def f1(x,y,i,j,k):
    return xdarr[i][j][k]

# second derivative of x
def f2(x,y,i,j,k):
    F = 0
    if (pointlist[j][k].left == False and pointlist[j][k].right == False):
        if (pointlist[j][k].top == True):
            F = - k/m * (xarr[i][j][k]-xarr[i][j][k-1]) + k/m * (xarr[i][j][k+1]-xarr[i][j][k]) \
                    + k/m * ((xarr[i][j+1][k+1]-xarr[i][j][k]) + (xarr[i][j+1][k-1] - xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
        elif (pointlist[j][k].bottom == True):
            F = - k/m * (xarr[i][j][k]-xarr[i][j][k-1]) + k/m * (xarr[i][j][k+1]-xarr[i][j][k]) \
                    + k/m * ((xarr[i][j-1][k+1]-xarr[i][j][k]) + (xarr[i][j-1][k-1]-xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
        else:
            F = - k/m * (xarr[i][j][k]-xarr[i][j][k-1]) + k/m * (xarr[i][j][k+1]-xarr[i][j][k]) \
                    + k/m * ((xarr[i][j+1][k+1] - xarr[i][j][k]) + (xarr[i][j+1][k-1] - xarr[i][j][k]) \
                    + (xarr[i][j-1][k+1]-xarr[i][j][k]) + (xarr[i][j-1][k-1]-xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]

    if (pointlist[j][k].left == True):
        if (pointlist[j][k].top == True):
            F = k/m * (xarr[i][j][k+1]-xarr[i][j][k]) + k/m * ((xarr[i][j+1][k+1] - xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
        elif (pointlist[j][k].bottom == True):
            F = k/m * (xarr[i][j][k+1]-xarr[i][j][k]) + k/m * ((xarr[i][j-1][k+1]-xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
        else:
            F = k/m * (xarr[i][j][k+1]-xarr[i][j][k]) + k/m * ((xarr[i][j+1][k+1] - xarr[i][j][k])  \
                    + (xarr[i][j-1][k+1]-xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]

    if (pointlist[j][k].right == True):
        if (pointlist[j][k].top == True):
            F = - k/m * (xarr[i][j][k]-xarr[i][j][k-1]) + k/m * ((xarr[i][j+1][k-1] - xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
        elif (pointlist[j][k].bottom == True):
            F = - k/m * (xarr[i][j][k]-xarr[i][j][k-1]) + k/m * ((xarr[i][j-1][k-1]-xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
        else:
            F = - k/m * (xarr[i][j][k]-xarr[i][j][k-1]) \
                    + k/m * ((xarr[i][j+1][k-1] - xarr[i][j][k]) + (xarr[i][j-1][k-1]-xarr[i][j][k])) \
                    - b/m * xdarr[i][j][k]
    return F

# first derivative of y
def f3(x,y,i,j,k):
    return ydarr[i][j][k]

# second derivative of y
def f4(x,y,i,j,k):
    F = 0
    if (pointlist[j][k].top == False and pointlist[j][k].bottom == False):
        if (pointlist[j][k].left == True):
            F = - k/m * (yarr[i][j][k]-yarr[i][j-1][k]) + k/m * (yarr[i][j+1][k]-yarr[i][j][k])\
                    + k/m * ((yarr[i][j+1][k+1] - yarr[i][j][k]) + (yarr[i][j-1][k+1]-yarr[i][j][k])) \
                    - b/m * ydarr[i][j][k]
        elif (pointlist[j][k].right == True):
            F = - k/m * (yarr[i][j][k]-yarr[i][j-1][k]) + k/m * (yarr[i][j+1][k]-yarr[i][j][k])\
                    + k/m * ((yarr[i][j+1][k-1] - yarr[i][j][k]) + (yarr[i][j-1][k-1]-yarr[i][j][k])) \
                    - b/m * ydarr[i][j][k]
        else:
            F = - k/m * (yarr[i][j][k]-yarr[i][j-1][k]) + k/m * (yarr[i][j+1][k]-yarr[i][j][k])\
                    + k/m * ((yarr[i][j+1][k+1] - yarr[i][j][k]) + (yarr[i][j+1][k-1] - yarr[i][j][k]) \
                    + (yarr[i][j-1][k+1]-yarr[i][j][k]) + (yarr[i][j-1][k-1]-yarr[i][j][k])) \
                    - b/m * ydarr[i][j][k]

    if (pointlist[j][k].top == True):
        if (pointlist[j][k].left == True):
            F = k/m * (yarr[i][j+1][k]-yarr[i][j][k])\
                    + k/m * ((yarr[i][j+1][k+1] - yarr[i][j][k])) - b/m * ydarr[i][j][k]
        elif (pointlist[j][k].right == True):
            F = k/m * (yarr[i][j+1][k]-yarr[i][j][k])\
                    + k/m * ((yarr[i][j+1][k-1] - yarr[i][j][k])) - b/m * ydarr[i][j][k]
        else:
            F = k/m * (yarr[i][j+1][k]-yarr[i][j][k])\
                    + k/m * ((yarr[i][j+1][k+1] - yarr[i][j][k]) + (yarr[i][j+1][k-1] - yarr[i][j][k])) \
                    - b/m * ydarr[i][j][k]

    if (pointlist[j][k].bottom == True):
        if (pointlist[j][k].left == True):
            F = - k/m * (yarr[i][j][k]-yarr[i][j-1][k]) + k/m * ((yarr[i][j-1][k+1]-yarr[i][j][k])) \
                    - b/m * ydarr[i][j][k]
        elif (pointlist[j][k].right == True):
            F = - k/m * (yarr[i][j][k]-yarr[i][j-1][k]) + k/m * (yarr[i][j-1][k-1]-yarr[i][j][k]) \
                    - b/m * ydarr[i][j][k]
        else:
            F = - k/m * (yarr[i][j][k]-yarr[i][j-1][k]) \
                    + k/m * ((yarr[i][j-1][k+1]-yarr[i][j][k]) + (yarr[i][j-1][k-1]-yarr[i][j][k])) \
                    - b/m * ydarr[i][j][k]

    return F

def position(steps, finalt, N, dt):
    global xarr
    global xdarr
    global yarr
    global ydarr

    tarr = np.linspace(0,finalt,steps)
    xarr = np.zeros(steps * N * N).reshape(steps,N,N)
    yarr = np.zeros(steps * N * N).reshape(steps,N,N)
    xdarr = np.zeros(steps * N * N).reshape(steps,N,N)
    ydarr = np.zeros(steps * N * N).reshape(steps,N,N)

    ydarr[0][N-1][N-1] = .5
    xdarr[0][N-1][N-1] = .5


    # gets initial conditions from each clothPoint
    for i in range(0,steps):
        for j in range(0,N):
            for k in range(0,N):
                xarr[0][j][k] = pointlist[j][k].xpos[0]
                yarr[0][j][k] = pointlist[j][k].ypos[0]

    for i in range(0,steps-1):
        for j in range(0,N):
            for k in range(0,N):
                xarr[i+1][j][k] = rk4(dt,tarr[i],xarr[i][j][k],f1,i,j,k)
                xdarr[i+1][j][k] = rk4(dt,tarr[i],xdarr[i][j][k],f2,i,j,k)
                yarr[i+1][j][k] = rk4(dt,tarr[i],yarr[i][j][k],f3,i,j,k)
                ydarr[i+1][j][k] = rk4(dt,tarr[i],ydarr[i][j][k],f4,i,j,k)

    # transpose position arrays
    x = np.transpose(xarr)
    y = np.transpose(yarr)

    # assigns the new position arrays to the clothPoint objects
    for i in range(0,steps-1):
        for j in range(0,N):
            for k in range(0,N):
                pointlist[j][k].xpos[i] = x[j][k][i]
                pointlist[j][k].ypos[i] = y[j][k][i]

    # stationary plots in 2D
    for j in range(0,N):
        for k in range(0,N):
            plt.plot(tarr,x[j][k])
            plt.xlabel('t')
            plt.ylabel('x')
    plt.show()
    for j in range(0,N):
        for k in range(0,N):
            plt.plot(tarr,y[j][k])
            plt.xlabel('t')
            plt.ylabel('y')
    plt.show()
    for j in range(0,N):
        for k in range(0,N):
            plt.plot(x[j][k],y[j][k])
            plt.xlabel('x')
            plt.ylabel('y')
    plt.show()
    
    return x, y

x, y = position(steps, finalt, N, dt)

""" Animation """
maxx = 0
maxy = 99
minx = 0
miny = 99

for j in range(0,N):
    for k in range(0,N):
        if j+max(x[j][k]) > maxx:
            maxx = j+max(x[j][k])
        if k+max(y[j][k]) > maxy:
            maxy = k+max(y[j][k])
        if j+min(x[j][k]) < minx:
            minx = j+min(x[j][k])
        if k+min(y[j][k]) < miny:
            miny = k+min(y[j][k])



fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', xlim=(minx-1, maxx+1),ylim=(miny-1, maxx+1))
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('test coupled oscillators')

lines = []

for i in range(N):
    line, = ax.plot([], [], 'o-')
    lines.append(line)

def init():
    for i in range(N):
        lines[i].set_data([], [])
    return lines

def animate(i):
    thisx = []
    thisy = []
    for j in range(N):
        for k in range(N):
            thisx.append([])
            thisy.append([])

            thisx[j].append(x[k,j][i]+1*k)
            thisy[j].append(y[k,j][i]+1*j)

    for j in range(N):

        lines[j].set_data(thisx[j],thisy[j])
    return lines

ani = animation.FuncAnimation(fig, animate, interval=10, blit=True, init_func=init)
plt.show()
