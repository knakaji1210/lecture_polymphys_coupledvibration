import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def doublePendulum(s, t):

    g = 9.80665             # [m/s^2] gravitational acceleration

    theta1, dtheta1, theta2, dtheta2 = s
    diff_theta = theta2 - theta1
    dtheta1dt = (m2*g*np.sin(theta2)*np.cos(diff_theta)
        + m2*l1*(dtheta1**2)*np.sin(diff_theta)*np.cos(diff_theta)
        + m2*l2*(dtheta2**2)*np.sin(diff_theta)
        - (m1 + m2)*g*np.sin(theta1))/(m1*l1 + m2*l1*(np.sin(diff_theta)**2))
    dtheta2dt = (-(m1 + m2)*l1*(dtheta1**2)*np.sin(diff_theta)
        - (m1 + m2)*g*np.sin(theta2)
        - m2*l2*(dtheta2**2)*np.sin(diff_theta)*np.cos(diff_theta)
        + (m1 + m2)*g*np.sin(theta1)*np.cos(diff_theta))/(m1*l2 + m2*l2*(np.sin(diff_theta)**2))
    return [dtheta1, dtheta1dt, dtheta2, dtheta2dt]

theta1_0 = np.pi        # [rad] initial angle
theta2_0 = np.pi/6      # [rad] initial angle
l1 = 1                  # [m] length of pendulum 
l2 = 1                  # [m] length of pendulum 
v1_0 = 0                # [m/s] initial velocity
v2_0 = 0                # [m/s] initial velocity
m1 = 1                  # [kg] mass
m2 = 1                  # [kg] mass
s0 = [theta1_0, v1_0/l1, theta2_0, v2_0/l2]     # initial condition
t = 10                  # [s] duration time
dt = 0.05               # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(doublePendulum, s0, t)
theta1, theta2 = sol[:,0], sol[:,2]
x1 = l1*np.sin(theta1)
y1 = -l1*np.cos(theta1)
x2 = x1 + l2*np.sin(theta2)
y2 = y1 - l2*np.cos(theta2)

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal', xlim=[-l1-l2, l1+l2], ylim=(-l1-l2, l1+l2))
ax.grid()

line, = ax.plot([], [], 'ro-', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line, time_text

def update(i):
    next_x = [0, x1[i], x2[i]]
    next_y = [0, y1[i], y2[i]]
    line.set_data(next_x, next_y)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/doublePendulum.gif', writer='pillow', fps=fps)