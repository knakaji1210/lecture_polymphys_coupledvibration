import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k, m):

    x, v_x = s              # s = (x, v_x)
    dsdt = [v_x, (-k*x)/m]
    return dsdt

k = 100                     # [N/m] spring constant
m = 50                      # [kg] mass
l = 20                      # [m] equilibrium length
x0 = 5                     # [m]   initial position
v_x0 = 0                    # [m/s] initial velocity
s0 = [x0, v_x0]             # initial condition
afreq = np.sqrt(k/m)        # angular frequency    
t = 4*np.pi/afreq           # [s] duration time
dt = 0.1                    # [s] interval time
ds = 0.01
smax = 6

amp = x0

t = np.arange(0, t+dt, dt)
s = np.arange(0, smax*np.pi, ds)

print(len(t))

sol = odeint(harmonicOscillator, s0, t, args=(k,m))
x = sol[:, 0]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.grid()

line, = plt.plot([], [], 'ro-', animated=True)

time_template = 'time = %.1fs'
time_text = fig.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
#    return line, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    ax.cla()
    ax.set_xlabel('x position [m]')
    ax.set_xlim(0,50)
    ax.set_ylim(-5,5)
    ax.grid()
    sp_x = [(l + amp*np.cos(4*np.pi*i/(len(t)-1)))*(np.sin((11/2)*s - np.pi/2)+(4/np.pi)*s + 1)/(smax*4+2) for s in s]
    sp_y = [0.5*np.cos((11/2)*s - np.pi/2) for s in s]
    ax.plot(sp_x,sp_y)
    ax.scatter(l + x[i], 0, c='red', s=200)
#    line.set_data([0, l + x[i]], [0, 0])
    time_text.set_text(time_template % (i*dt))
#    return line, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, interval=frame_int, repeat=True)

plt.show()

ani.save('./gif/harmonicOscillator_spring.gif', writer='pillow', fps=fps)