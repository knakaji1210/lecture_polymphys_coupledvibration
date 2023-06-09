# ordinary differential equation of damped harmonic oscillator

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedHarmonicOscillator(s, t, k, m, c):

    x, v = s                # s = (x, v)
    dsdt = [v, (-k*x)/m - c*v/m]    # 粘性抵抗の項を追加
    return dsdt

# variables
try:
    k = float(input('spring constant [N/m] (default=10.0): '))
except ValueError:
    k = 10.0                # [N/m] spring constant
try:
    m = float(input('mass [kg] (default=1.0): '))
except ValueError:
    m = 1.0                 # [kg] mass
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5                 # [m/s] damping coefficient
l = 20                      # [m] equilibrium length
afreq0 = np.sqrt(k/m)       # natural　angular frequency
rho = c/(2*m)
tau = 1/rho
if np.abs(rho - afreq0) < 0.05: # critical damping (rho == afreq0はほぼ無理)
    period = tau                # period [s] (T)
    cond = "cd"
elif rho - afreq0 < 0: # under damping
    afreq = np.sqrt(afreq0**2 - rho**2)
    period = 2*np.pi/afreq      # period [s] (T)
    cond = "ud"
elif rho - afreq0 > 0: # over damping
    eta = np.sqrt(rho**2 - afreq0**2)
    period = 1/(rho - eta)      # period [s] (T)
    cond = "od"

tmax = 6*period             # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
try:
    x0 = float(input('initial position [m] (default=5.0): '))
except ValueError:
    x0 = 5.0
try:
    v0 = float(input('initial velocity [m/s] (default=0.0): '))
except ValueError:
    v0 = 0.0

s0 = [x0, v0]               # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(dampedHarmonicOscillator, s0, t, args=(k,m,c))
x = sol[:, 0]

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-5, 50), ylim=(-5, 5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
var_template = r'$k$ = {0:.1f} N/m, $m$ = {1:.1f} kg, $c$ = {2:.1f} kg/s'.format(k,m,c)
ax.text(0.4, 0.9, var_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
period_template = r'$\tau$ = {0:.2f} s, $T$ = {1:.2f} s ({2})'.format(tau,period,cond)
ax.text(0.1, 0.8, period_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

line, = plt.plot([], [], 'ro-', animated=True)
# ここでは[],[]としているが、下でline.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return line, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
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

savefile = './gif/damped_harmonicOsci_(x={0:.1f},v={1:.1f},k={2:.1f},m={3:.1f},c={4:.1f}).gif'.format(x0,v0,k,m,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()