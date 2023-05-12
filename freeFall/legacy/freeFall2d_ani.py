import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def freeFall(s, t, r, m):

    g  = 9.80665            # [m/s^2] gravitational acceleration
    x, y, v_x, v_y = s  # s = (x, y, v_x, v_y)
    v = np.sqrt(v_x*v_x + v_y*v_y)
    dsdt = [v_x, v_y, 0, (-r*v-m*g)/m]
    return dsdt

r = 50                      # air resistance
m = 100                     # [kg] mass
x0 = 0                      # [m]   initial position of x
y0 = 0                      # [m]   initial position of y
v_x0 = 10                   # [m/s] initial velocity of x
v_y0 = 20                   # [m/s] initial velocity of y
s0 = [x0, y0, v_x0, v_y0]   # initial condition
t = 5                      # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(freeFall, s0, t, args=(r,m))
x, y = sol[:, 0], sol[:, 1]

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', ylabel='y position [m]', xlim=[-10, 100], ylim=(-200, 50))
ax.grid()

ball, = ax.plot([], [], 'ro', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return ball, time_text

def update(i):
    ball.set_data(x[i], y[i])
    time_text.set_text(time_template % (i*dt))
    return ball, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/freeFall2d.gif', writer='pillow', fps=fps)