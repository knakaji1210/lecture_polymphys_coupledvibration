# animation of sinusoidal vibration (spring)

import numpy as np
import matplotlib.pyplot as plt
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
    x_tri = [(l/4)+((l/4)*np.sin(afreq*t[i])+l/2)*s for s in s]
    y_tri = [(2/3)*np.arccos(np.cos(8*np.pi*s-np.pi/2+0.1))-1 for s in s]
    x_rod1 = [(l/4)*s for s in s]
    x_rod2 = [(3*l/4)+(l/4)*s+amp*np.sin(afreq*t[i]) for s in s]
    y_rod = np.zeros(len(s))
    x_mass = l + amp*np.sin(afreq*t[i])
    y_mass = 0
    ax.plot(x_tri,y_tri, c='b', zorder=1)
    ax.plot(x_rod1,y_rod, c='b', zorder=2)
    ax.plot(x_rod2,y_rod, c='b', zorder=3)
    ax.scatter(x_mass,y_mass, s=100, c='r', zorder=4)        
    
f = np.arange(0, len(t))
frame_int = 100 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, plot, frames=f, interval=frame_int, repeat=False)

ani.save('./gif/spring.gif', writer='pillow', fps=fps)

plt.show()