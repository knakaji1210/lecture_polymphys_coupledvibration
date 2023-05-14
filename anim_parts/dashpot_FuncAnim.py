# animation of sinusoidal vibration (dashpot parts)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

# variables
k = 10                      # [N/m] spring constant
m = 50                      # [kg] mass
l = 20                      # [m] equilibrium length
amp = 5                     # [m] vibration amplitude
afreq = np.sqrt(k/m)        # angular frequency
period = 2*np.pi/afreq      # period [s] (T)
tmax = 2*period             # [s] duration time (2*T)
dt = tmax/100               # [s] interval time
smax = 1                    # [m]
ds = smax/1000              # [m]

t = np.arange(0, tmax+dt, dt)
s = np.arange(0, smax, ds)

fig = plt.figure()
ax = fig.add_subplot(111)

def plot(i):
    ax.cla()
    ax.set_xlabel('$x$ position [m]')
    ax.set_xlim(-5,50)
    ax.set_ylim(-5,5)
    ax.grid()
    ax.set_axisbelow(True)
    x_mass = l + amp*np.sin(afreq*t[i])
    y_mass = 0
    x_d0 = [(l/4)*s for s in s]
    x_d1 = [(l/4)*(1+2.1*s) for s in s]
    x_d2 = [(l/2)+(l/2)*s+amp*np.sin(afreq*t[i]) for s in s]
    y_d0 = np.zeros(len(s))
    h = 0.4
    y_d1 = h*np.ones(len(s))
    y_d2 = -h*np.ones(len(s))
    rect = patches.Rectangle(xy=(l/4, -h), width=l/2, height=2*h, facecolor='y')
    ax.plot([l/4, l/4],[h,-h], c='b', zorder=1)
    ax.plot([(l/2)+amp*np.sin(afreq*t[i]), (l/2)+amp*np.sin(afreq*t[i])],[0.7*h,-0.7*h], lw=4, c='b', zorder=2)
    ax.plot(x_d0,y_d0, c='b', zorder=3)
    ax.plot(x_d1,y_d1, c='b', zorder=4)
    ax.plot(x_d1,y_d2, c='b', zorder=5)
    ax.plot(x_d2,y_d0, c='b', zorder=6)
    ax.add_patch(rect)
    ax.scatter(x_mass,y_mass, s=100, c='r', zorder=7)     

f = np.arange(0, len(t))
frame_int = 100 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, plot, frames=f, interval=frame_int, repeat=False)

ani.save('./gif/dashpot.gif', writer='pillow', fps=fps)

plt.show()