import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k, m):

    x, v = s                # s = (x, v)
    dsdt = [v, (-k*x)/m - c*v/m]
    return dsdt

k = 100                     # [N/m] spring constant
m = 10                      # [kg] mass
c = 10                      # air resistance
l = 20                      # [m] equilibrium length
x0 = 5                      # [m]   initial position
v0 = 0                    # [m/s] initial velocity
s0 = [x0, v0]             # initial condition
afreq = np.sqrt(k/m)        # angular frequency    
t = 6*np.pi/afreq           # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(harmonicOscillator, s0, t, args=(k,m))
x = sol[:, 0]

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', ylabel='y position [m]', xlim=(0, 40), ylim=(-10, 10))
ax.grid()

line, = plt.plot([], [], 'ro-', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l + x[i]], [0, 0])
    time_text.set_text(time_template % (i*dt))
    return line, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/dampled_harmonicOscillator.gif', writer='pillow', fps=fps)