# animation of sinusoidal vibration (dashpot)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
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
w = 0.4                     # ratio of dashpot width

t = np.arange(0, tmax+dt, dt)
s = np.arange(0, smax, ds)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=[-5,50], ylim=(-5,5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
x_d0 = [(l/4)*s for s in s]
y_d0 = np.zeros(len(s))
x_d1 = [(l/4)*(1+2.1*s) for s in s]
y_d1 = w*np.ones(len(s))
rect = patches.Rectangle(xy=(l/4, -w), width=l/2, height=2*w, facecolor='y')
ax.plot([l/4, l/4],[w,-w], c='b')
ax.plot(x_d0,y_d0, c='b')
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,-y_d1, c='b')
ax.add_patch(rect)

rod, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
mass, = ax.plot([],[], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_dataで実際の値を入れている

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():
    time_text.set_text('')
    return rod, damper, mass, time_text

def update(i):
    x_rod = [(l/2)+(l/2)*s+0.9*amp*np.sin(afreq*t[i]) for s in s]
    x_damp = (l/2)+0.9*amp*np.sin(afreq*t[i])
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod.set_data(x_rod,y_d0)
    x_mass = [0, l + amp*np.sin(afreq*t[i])]
    y_mass = [0, 0]
    damper.set_data(x_damper,y_damper)
    time_text.set_text(time_template % (i*dt))
    mass.set_data(x_mass,y_mass)
    return rod, damper, mass, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

ani.save('./gif/dashpot.gif', writer='pillow', fps=fps)

plt.show()