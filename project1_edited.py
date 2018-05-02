'''
Alexander Bart and Julie Bosley
Feb 19 2018
PHGN498A - Project 1

N x N Coupled Pendula
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import mpl_toolkits.mplot3d.axes3d as p3
'''
Coupled Pendula Angles
'''
def coupled_pendula_3D(steps,timef,N):

    finaltime = timef
    dt = finaltime/steps
    theta = np.zeros((N,N,steps+1))
    phi = np.zeros((N,N,steps+1))

    thdot = np.zeros((N,N,steps+1))
    phdot = np.zeros((N,N,steps+1))

    thdot[1,1,0] = 2
    phdot[1,1,0] = -1

    k = 10
    m = 1
    R = 5
    g = 98.1

    for i in range(steps):
        for j in range(N):
            if j == 0:
                thdot[:,j,i+1] = thdot[:,j,i] + dt*(-g/R*np.sin(theta[:,j,i]) +\
                k/m*np.cos(theta[:,j,i])*(np.sin(theta[:,j+1,i])-np.sin(theta[:,j,i])))

                phdot[j,:,i+1] = phdot[j,:,i] + dt*(-g/R*np.sin(phi[j,:,i]) +\
                k/m*np.cos(phi[j,:,i])*(np.sin(phi[j+1,:,i])-np.sin(phi[j,:,i])))

            if j == N-1:
                thdot[:,j,i+1] = thdot[:,j,i] + dt*(-g/R*np.sin(theta[:,j,i]) +\
                k/m*np.cos(theta[:,j,i])*(np.sin(theta[:,j-1,i])-np.sin(theta[:,j,i])))

                phdot[j,:,i+1] = phdot[j,:,i] + dt*(-g/R*np.sin(phi[j,:,i]) +\
                k/m*np.cos(phi[j,:,i])*(np.sin(phi[j-1,:,i])-np.sin(phi[j,:,i])))

            else:
                thdot[:,j,i+1] = thdot[:,j,i] + dt*(-g/R*np.sin(theta[:,j,i]) +\
                k/m*np.cos(theta[:,j,i])*\
                (np.sin(theta[:,j+1,i])-2*np.sin(theta[:,j,i])+theta[:,j-1,i]))

                phdot[j,:,i+1] = phdot[j,:,i] + dt*(-g/R*np.sin(phi[j,:,i]) +\
                k/m*np.cos(phi[j,:,i])*\
                (np.sin(phi[j+1,:,i])-2*np.sin(phi[j,:,i])+phi[j-1,:,i]))

            theta[:,j,i+1] = theta[:,j,i] + dt*thdot[:,j,i+1]
            phi[j,:,i+1] = phi[j,:,i] + dt*phdot[j,:,i+1]

    return [theta,phi]


steps = 10000
finaltime = 10
N= 3
dt = finaltime/steps
time = np.linspace(0, finaltime, steps + 1)

R = 5

angles = coupled_pendula_3D(steps,finaltime,N)

x = np.zeros((N,N,steps+1))
y = np.zeros((N,N,steps+1))
z = np.zeros((N,N,steps+1))

for i in range(N):
    for j in range(N):
        x[i,j] = R * np.sin(angles[0][i,j]) + j*R
        y[i,j] = R * np.sin(angles[1][i,j]) + i*R
        z[i,j] = -R * np.cos(angles[1][i,j])*np.cos(angles[0][i,j])

fig = plt.figure()
ax = p3.Axes3D(fig)

ax.set_xlim3d([-R-1, N*R+1])
ax.set_xlabel('X')

ax.set_ylim3d([-R-1, N*R+1])
ax.set_ylabel('Y')

ax.set_zlim3d([-R-1, 1])
ax.set_zlabel('Z')

ax.set_title('3D NxN Coupled Pendulums')

lines = []

for i in range(N):
    line, = ax.plot([], [], [], 'o')
    lines.append(line)

def init():
    for i in range(N):
        lines[i].set_data([], [])
        lines[i].set_3d_properties([])

    return lines

def animate(i):
    thisx = []
    thisy = []
    thisz = []
    for j in range(N):
        for k in range(N):
            thisx.append([])
            thisy.append([])
            thisz.append([])
            
            thisx[j].append(x[k,j][i])
            thisy[j].append(y[k,j][i])
            thisz[j].append(z[k,j][i])

    for j in range(N):
        #set_data() only sets 2 arrays
        lines[j].set_data(thisx[j],thisy[j])
        #need this command to set third dat set
        lines[j].set_3d_properties(thisz[j])

    return lines

ani = animation.FuncAnimation(fig, animate, interval=1, blit=True, init_func=init)
plt.show()

