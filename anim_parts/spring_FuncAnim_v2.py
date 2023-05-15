# animation of sinusoidal vibration (spring)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# variables
k = 100                     # [N/m] spring constant
m = 20                      # [kg] mass
l = 20                      # [m] equilibrium length
amp = 5                     # [m] vibration amplitude
afreq = np.sqrt(k/m)        # angular frequency
period = 2*np.pi/afreq      # period [s] (T)
tmax = 2*period             # [s] duration time (2*T)
dt = 0.05                   # [s] interval time
smax = 1                    # [m]
ds = smax/1000              # [m]

t = np.arange(0, tmax, dt)
s = np.arange(0, smax, ds)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=[-5,50], ylim=(-5,5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
x = [(l/4)*s for s in s]
y = np.zeros(len(s))
ax.plot(x,y, c='b')

rod, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
mass, = ax.plot([],[], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():
    time_text.set_text('')
    return rod, triangle, mass, time_text

def update(i):
    x_rod = [(3*l/4)+(l/4)*s+amp*np.sin(afreq*t[i]) for s in s]
    x_tri = [(l/4)+((l/4)*np.sin(afreq*t[i])+l/2)*s for s in s]
    y_tri = [(2/3)*np.arccos(np.cos(8*np.pi*s-np.pi/2+0.1))-1 for s in s]
    x_mass = [0, l + amp*np.sin(afreq*t[i])]
    y_mass = [0, 0]
    rod.set_data(x_rod,y)
    triangle.set_data(x_tri,y_tri)
    mass.set_data(x_mass,y_mass)
    time_text.set_text(time_template % (i*dt))
    return rod, triangle, mass, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

ani.save('./gif/spring.gif', writer='pillow', fps=fps)

plt.show()