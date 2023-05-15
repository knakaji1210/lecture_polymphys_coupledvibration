# animation of sinusoidal vibration (rod)

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
    x1 = [(l/4)*s for s in s]
    x2 = [(3*l/4)+(l/4)*s+amp*np.sin(afreq*t[i]) for s in s]
    y = np.zeros(len(s))
    ax.plot(x1,y, c='b')
    ax.plot(x2,y, c='b')
    

f = np.arange(0, len(t))
frame_int = 100 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, plot, frames=f, interval=frame_int, repeat=False)

ani.save('./gif/rod.gif', writer='pillow', fps=fps)

plt.show()