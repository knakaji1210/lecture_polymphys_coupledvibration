import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def freeFall(s, t, r, m):

    g  = 9.80665            # [m/s^2] gravitational acceleration
    y, v_y = s              # s = (y, v_y)
    dsdt = [v_y, (-r*v_y-m*g)/m]
    return dsdt

r = 50                      # air resistance
m = 100                     # [kg] mass
y0 = 0                      # [m]   initial position
v_y0 = 20                     # [m/s] initial velocity
s0 = [y0, v_y0]               # initial condition
t = 7                       # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(freeFall, s0, t, args=(r,m))
y = sol[:, 0]

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', ylabel='y position [m]', xlim=(-10, 10), ylim=(-100, 50))
ax.grid()

ball, = plt.plot([], [], 'ro', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return ball, time_text

def update(i):
    ball.set_data(0, y[i])
    time_text.set_text(time_template % (i*dt))
    return ball, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/freeFall.gif', writer='pillow', fps=fps)