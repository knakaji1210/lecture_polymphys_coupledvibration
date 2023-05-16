# animation of sinusoidal vibration (dashpot parts)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.animation import FuncAnimation

# variables
k = 100                     # [N/m] spring constant
m = 20                      # [kg] mass
l = 20                      # [m] equilibrium length
afreq = np.sqrt(k/m)        # angular frequency
period = 2*np.pi/afreq      # period [s] (T)
tmax = 2*period             # [s] duration time (2*T)
dt = 0.05                   # [s] interval time
w = 0.4                     # ratio of dashpot width

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
x_rod1 = [0, l/4]
y = [0, 0]
x_d1 = [l/4, 3.1*l/4]
y_d1 = [w, w]
y_d2 = [-w, -w]
ax.plot(x_rod1,y, c='b')
ax.plot(x_d1,y_d1, c='b')
ax.plot(x_d1,y_d2, c='b')
ax.plot([l/4,l/4],[w,-w], c='b')
rect = patches.Rectangle(xy=(l/4, -w), width=l/2, height=2*w, facecolor='y')
ax.add_patch(rect)

rod, = ax.plot([],[], 'b', animated=True)
damper, = ax.plot([],[], 'b', lw=4, animated=True)
# ここでは[],[]としているが、下でrod.set_dataで実際の値を入れている

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():
    time_text.set_text('')
    return rod, damper, time_text

def update(i):
    x_rod2 = [l/2 + 0.9*x[i], l + 0.9*x[i]]
    x_damp = (l/2)+0.9*x[i]
    y_damp = 0.7*w
    x_damper = [x_damp, x_damp]
    y_damper = [y_damp, -y_damp]
    rod.set_data(x_rod2,y)
    damper.set_data(x_damper,y_damper)
    time_text.set_text(time_template % (i*dt))
    return rod, damper, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f, 
                    init_func=init, blit=True, interval=frame_int, repeat=False)

ani.save('./gif/dashparts.gif', writer='pillow', fps=fps)

plt.show()