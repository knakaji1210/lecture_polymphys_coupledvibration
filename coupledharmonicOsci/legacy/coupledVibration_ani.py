import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k1, k2, m):

    x1, v1, x2, v2 = s
    a1 = - (k1+k2)*x1/m + k2*x2/m
    a2 = k2*x1/m - (k1+k2)*x2/m
    dsdt = [v1, a1, v2, a2]
    return dsdt

k1 = 10                    # [N/m] spring constant
k2 = 50                     # [N/m] spring constant
m = 10                      # [kg] mass
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
L = l1 + l2 + l3
x1_0 = 5                    # [m]   initial position
x2_0 = 10                    # [m]   initial position
v1_0 = 0                    # [m/s] initial velocity
v2_0 = 0                    # [m/s] initial velocity
s0 = [x1_0, v1_0, x2_0, v2_0]
af1 = np.sqrt(k1/m)         # angular frequency   
af2 = np.sqrt((k1+2*k2)/m)  # angular frequency   
t = 8*np.pi/af1             # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(harmonicOscillator, s0, t, args=(k1, k2, m))
x1, x2 = sol[:, 0], sol[:, 2]
q1 = (1/np.sqrt(2))*(x1 + x2)
q2 = (1/np.sqrt(2))*(-x1 + x2)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', xlim=(0, L), ylim=(-1, 3))
ax.grid()

line1, = plt.plot([], [], 'ro-', animated=True)
line2, = plt.plot([], [], 'bo-', animated=True)
line3, = plt.plot([], [], 'go-', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
#    return line1, line2, line3, time_text
    return line1, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line1.set_data([0, l1 + x1[i], l1 + l2 + x2[i], L], [0, 0, 0, 0])
#    line2.set_data([0, L/2 + q1[i]], [1, 1])
#    line3.set_data([0, L/2 + q2[i]], [2, 2])
    time_text.set_text(time_template % (i*dt))
#    return line1, line2, line3, time_text
    return line1, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/coupledVibration.gif', writer='pillow', fps=fps)