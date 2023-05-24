# ordinary differential equation of harmonic oscillator

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k, m):

    x, v = s                # s = (x, v)
    dsdt = [v, (-k*x)/m]    # dx/dt = v, m*dv/dt = -k*x 
    return dsdt

# variables
k = 100                     # [N/m] spring constant
m = 20                      # [kg] mass
l = 20                      # [m] equilibrium length
afreq = np.sqrt(k/m)        # angular frequency    
period = 2*np.pi/afreq      # period [s] (T)
t = 2*period                # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
try:
    x0 = int(input('initial position (default=5): '))
except ValueError:
    x0 = 5
try:
    v0 = int(input('initial velocity (default=0): '))
except ValueError:
    v0 = 0

s0 = [x0, v0]               # initial condition

t = np.arange(0, t+dt, dt)

sol = odeint(harmonicOscillator, s0, t, args=(k,m))  # ODEの解を求めている
#print(sol.shape) # (len(t), 2)が出てくる。sol = [[x],[v]]ということ
x = sol[:, 0]    # [x]が出てくる
#print(x.shape)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-5, 50), ylim=(-5, 5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
# ここでは[],[]としているが、下でline.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return line, time_text

def update(i):              # ここのiは下のframes=fに対応した引数になっている
    x_pos = [0, l + x[i]]
    y_pos = [0, 0]
    line.set_data(x_pos,y_pos)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/harmonicOsci_(x={0},v={1}).gif'.format(x0,v0)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()