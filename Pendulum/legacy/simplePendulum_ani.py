import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def simplePendulum(s, t):

    g = 9.80665             # [m/s^2] gravitational acceleration

    theta, dtheta = s
    dsdt = [dtheta, -(g/l)*np.sin(theta)]
    return dsdt

theta0 = np.pi/4                # [m]   initial angle
l = 1                   # [m] length of pendulum 
v0 = 1                    # [m/s] initial velocity
s0 = [theta0, v0/l]             # initial condition
t = 10                      # [s] duration time
dt = 0.05                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(simplePendulum, s0, t)
theta = sol[:,0]
x = l*np.sin(theta)
y = -l*np.cos(theta)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', xlim=[-l, l], ylim=(-l, l))
ax.grid()

line, = ax.plot([], [], 'ro-', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line, time_text

def update(i):
    next_x = [0, x[i]]
    next_y = [0, y[i]]
    line.set_data(next_x, next_y)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/simplePendulum.gif', writer='pillow', fps=fps)