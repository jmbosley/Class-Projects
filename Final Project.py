import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import random
import time
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axes3d as p3
import scipy.integrate as scint
import matplotlib.patches as patches

###############################################################################
#all with positive and negative
#competition coefficient
#represents the effect that one species has on the other

alpha12 = -2
alpha13 = 5
alpha14 = 0

alpha21 = 1
alpha23 = -2
alpha24 = 0

alpha31 = -2
alpha32 = 10
alpha34 = 0

alpha41 = 0
alpha42 = 0
alpha43 = 0
###############################################################################
##All equally competeing
#alpha12 = 1
#alpha13 = 1
#alpha14 = 0
#
#alpha21 = 1
#alpha23 = 1
#alpha24 = 0
#
#alpha31 = 1
#alpha32 = 1
#alpha34 = 0
#
#alpha41 = 0
#alpha42 = 0
#alpha43 = 0
###############################################################################
###pop 1 grows first then pop 2 then pop 3
#alpha12 = 0
#alpha13 = 0
#alpha14 = 0
#
#alpha21 = 1
#alpha23 = 0
#alpha24 = 0
#
#alpha31 = 1
#alpha32 = 1
#alpha34 = 0
#
#alpha41 = 0
#alpha42 = 0
#alpha43 = 0
###############################################################################

#intraspecific competition with population models
# compeition between individuals in the same species
steps = 10**2
tfinal = 10

#number of individuals
N1 = np.zeros(steps)
N1[0] = 10
N2 = np.zeros(steps)
N2[0] = 10
N3 = np.zeros(steps)
N3[0] = 10
N4 = np.zeros(steps)
N4[0] = 10

#time stuff
t = np.linspace(0,tfinal,steps)
dt = tfinal/steps

#rate of increase
r1 = 10
r2 = 10
r3 = 10
r4 = 10



#Carrying Capacity
k1 = 100
k2 = 100
k3 = 100
k4 = 100

#Runge Kutta
##################################################################################
#subscripts are flipped
def f1(N1,N2,N3,N4,t):
    return r1*N1*(k1-N1-alpha12*N2-alpha13*N3-alpha14*N4)/k1
def f2(N1,N2,N3,N4,t):
    return r2*N2*(k2-N2-alpha21*N1-alpha23*N3-alpha24*N4)/k2
def f3(N1,N2,N3,N4,t):
    return r3*N3*(k3-N3-alpha31*N1-alpha32*N2-alpha34*N4)/k3
def f4(N1,N2,N3,N4,t):
    return r4*N4*(k4-N4-alpha41*N1-alpha42*N2-alpha43*N3)/k4


def rk4(f1, f2, f3, t):
    
    for i in range(steps-1):
        k11 = dt*f1(N1[i],N2[i],N3[i],N4[i],t[i])
        k12 = dt*f2(N1[i],N2[i],N3[i],N4[i],t[i])
        k13 = dt*f3(N1[i],N2[i],N3[i],N4[i],t[i])
        k14 = dt*f3(N1[i],N2[i],N3[i],N4[i],t[i])
        
        k21 = dt*f1(N1[i]+k11/2,N2[i]+k12/2,N3[i]+k13/2,N4[i]+k14/2,t[i]+dt/2)
        k22 = dt*f2(N1[i]+k11/2,N2[i]+k12/2,N3[i]+k13/2,N4[i]+k14/2,t[i]+dt/2)
        k23 = dt*f3(N1[i]+k11/2,N2[i]+k12/2,N3[i]+k13/2,N4[i]+k14/2,t[i]+dt/2)
        k24 = dt*f4(N1[i]+k11/2,N2[i]+k12/2,N3[i]+k13/2,N4[i]+k14/2,t[i]+dt/2)
        
        k31 = dt*f1(N1[i]+k21/2,N2[i]+k22/2,N3[i]+k23/2,N4[i]+k24/2,t[i]+dt/2)
        k32 = dt*f2(N1[i]+k21/2,N2[i]+k22/2,N3[i]+k23/2,N4[i]+k24/2,t[i]+dt/2)
        k33 = dt*f3(N1[i]+k21/2,N2[i]+k22/2,N3[i]+k23/2,N4[i]+k24/2,t[i]+dt/2)
        k34 = dt*f4(N1[i]+k21/2,N2[i]+k22/2,N3[i]+k23/2,N4[i]+k24/2,t[i]+dt/2)
        
        k41 = dt*f1(N1[i]+k31/2, N2[i]+k32/2,N3[i]+k33/2,N4[i]+k34/2, t[i]+dt/2)
        k42 = dt*f2(N1[i]+k31/2, N2[i]+k32/2,N3[i]+k33/2,N4[i]+k34/2, t[i]+dt/2)
        k43 = dt*f3(N1[i]+k31/2, N2[i]+k32/2,N3[i]+k33/2,N4[i]+k34/2, t[i]+dt/2)
        k44 = dt*f4(N1[i]+k31/2, N2[i]+k32/2,N3[i]+k33/2,N4[i]+k34/2, t[i]+dt/2)
        
        N1[i+1] = N1[i] + dt/6. * (k11 + 2. * k21 + 2. * k31 + k41)
        N2[i+1] = N2[i] + dt/6. * (k12 + 2. * k22 + 2. * k32 + k42)
        N3[i+1] = N3[i] + dt/6. * (k13 + 2. * k23 + 2. * k33 + k43)
        N4[i+1] = N4[i] + dt/6. * (k14 + 2. * k24 + 2. * k34 + k44)
    return (N1,N2,N3,N4)


N1,N2,N3,N4 = rk4(f1,f2,f3,t)
#animation starts here
##################################################################################
fig = plt.figure()
ax = fig.add_subplot(121, aspect='auto', xlim=(-1, tfinal),ylim=(-5, k1+5))
ax.set_xlabel('Time')
ax.set_ylabel('Population')
ax.set_title('Population vs Time')

line1, = ax.plot([], [], '-or')
line2, = ax.plot([], [], '-ob')
line3, = ax.plot([], [], '-og')
line4, = ax.plot([], [], '-ok')
ax.legend(handles=[line4, line1,line2,line3], labels=['Population 0','Population 1','Population 2','Population 3'],loc=2)

 
def animate(i):
    line1.set_data(t[:i], N1[:i])
    line2.set_data(t[:i], N2[:i])
    line3.set_data(t[:i], N3[:i])
    line4.set_data(t[:i], N4[:i])
    
    return line1, line2, line3, line4
     
anim = animation.FuncAnimation(fig, animate,  interval=500, repeat=True,blit = True)
##################################################################################
ax1 = fig.add_subplot(122, aspect='auto', xlim=(-100, 100), ylim=(-100, 100))
ax1.set_yticklabels([])
ax1.set_xticklabels([])

scatter, = ax1.plot([],[], 's', color='r')
scatter2, = ax1.plot([],[],'s',color = 'b')
scatter3, = ax1.plot([], [], 's', color = 'g')

def animate2(i):
    xs = []
    ys = []
    xs2 = []
    ys2 = []
    xs3 = []
    ys3 = []
    maxnum = int(N1[i])
    maxnum2 = int(N2[i])
    maxnum3 = int(N3[i])
    for k in range(maxnum):
        xpos = float(random.randint(-100, 100))
        ypos = float(random.randint(-100, 100))
        xs.append(xpos)
        ys.append(ypos)
    for j in range(maxnum2):
        xpos2 = float(random.randint(-100, 100))
        ypos2 = float(random.randint(-100, 100))
        xs2.append(xpos2)
        ys2.append(ypos2)
    for j in range(maxnum3):
        xpos3 = float(random.randint(-100, 100))
        ypos3 = float(random.randint(-100, 100))
        xs3.append(xpos3)
        ys3.append(ypos3)     
    scatter.set_data(xs, ys)
    scatter2.set_data(xs2, ys2)
    scatter3.set_data(xs3, ys3)  
    return scatter, scatter2, scatter3, #scatter4
anim2 = animation.FuncAnimation(fig, animate2, interval=500, blit=True, repeat=True)
##################################################################################
plt.show()



