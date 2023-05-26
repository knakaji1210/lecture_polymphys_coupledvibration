import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k, m):

    x1, v1, x2, v2, x3, v3 = s
    a1 = -2*k*x1/m + k*x2/m  - c*v1/m
    a2 = k*x1/m - 2*k*x2/m + k*x3/m  - c*v2/m
    a3 = k*x2/m - 2*k*x3/m  - c*v3/m
    dsdt = [v1, a1, v2, a2, v3, a3]
    return dsdt

k = 100                      # [N/m] spring constant
m = 10                      # [kg] mass
c = 5                      # air resistance
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
l4 = 20                     # [m] equilibrium length
L = l1 + l2 + l3 + l4
x1_0 = 10                    # [m]   initial position
x2_0 = 0                   # [m]   initial position
x3_0 = 0                   # [m]   initial position
v1_0 = 0                    # [m/s] initial velocity
v2_0 = 0                    # [m/s] initial velocity
v3_0 = 0                    # [m/s] initial velocity
s0 = [x1_0, v1_0, x2_0, v2_0, x3_0, v3_0]
af = np.sqrt(2*k/m)         # angular frequency    
t = 16*np.pi/af             # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(harmonicOscillator, s0, t, args=(k, m))
x1, x2, x3 = sol[:, 0], sol[:, 2], sol[:, 4]
q1 = (1/2)*(x1 + np.sqrt(2)*x2+x3)
q2 = (1/2)*(-np.sqrt(2)*x1 + np.sqrt(2)*x3)
q3 = (1/2)*(x1 - np.sqrt(2)*x2+x3)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', xlim=(0, L), ylim=(-L/4, L/2))
ax.grid()

line1, = plt.plot([], [], 'ro-', animated=True)
line2, = plt.plot([], [], 'bo-', animated=True)
line3, = plt.plot([], [], 'go-', animated=True)
line4, = plt.plot([], [], 'yo-', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line1, line2, line3, line4, time_text
#    return line1, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line1.set_data([0, l1, l1 + l2, l1 + l2 + l3, L], [0, x1[i], x2[i], x3[i], 0])
    line2.set_data([0, L/2 + q1[i]], [10, 10])
    line3.set_data([0, L/2 + q2[i]], [20, 20])
    line4.set_data([0, L/2 + q3[i]], [30, 30])
    time_text.set_text(time_template % (i*dt))
    return line1, line2, line3, line4, time_text
#    return line1, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/damped_tripledVibration.gif', writer='pillow', fps=fps)