# animation of sinusoidal vibration (triagnle wave)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# variables
k = 100                     # [N/m] spring constant
m = 20                      # [kg] mass
l = 20                      # [m] equilibrium length
afreq = np.sqrt(k/m)        # angular frequency
period = 2*np.pi/afreq      # period [s] (T)
tmax = 2*period             # [s] duration time (2*T)
dt = 0.05                   # [s] interval time

# initial condition
try:
    x0 = int(input('initial position (default=5): '))
except ValueError:
    x0 = 5                  # [m] vibration amplitude

t = np.arange(0, tmax, dt)

x = [x0*np.cos(afreq*tim) for tim in t]

fig = plt.figure()
ax = fig.add_subplot(111, xlim=[-5,50], ylim=(-5,5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
y = [0, 0]

triangle, = ax.plot([],[], 'b', animated=True)
# ここでは[],[]としているが、下でtriangle.set_dataで実際の値を入れている

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():
    time_text.set_text('')
    return triangle, time_text

def update(i):
    x_tri = np.linspace(l/4, 3*l/4 + x[i],100)
    y_tri = (2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(x[i]+l/2)-np.pi/2+0.1))-1
    triangle.set_data(x_tri,y_tri)
    time_text.set_text(time_template % (i*dt))
    return triangle, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

ani.save('./gif/triangle.gif', writer='pillow', fps=fps)

plt.show()