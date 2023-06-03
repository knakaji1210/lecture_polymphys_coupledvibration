# ordinary differential equation of damped harmonic oscillator

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def forcedHarmonicOscillator(s, t, k, m, c, Fa, af):

    x, v = s                # s = (x, v)
    dsdt = [v, (-k*x)/m - c*v/m + (Fa/m)*np.sin(af*t)]    # 粘性抵抗、強制振動の項を追加
    return dsdt

# variables
try:
    k = float(input('spring constant [N/m] (default=10.0): '))
except ValueError:
    k = 10.0               # [N/m] spring constant
try:
    m = float(input('mass [kg] (default=1.0): '))
except ValueError:
    m = 1.0                # [kg] mass
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5
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

tmax = 10*period             # [s] duration time
dt = 0.05                    # [s] interval time

# forced oscillation
try:
    Fa = float(input('amplitude for forced oscillation [N] (default=10.0): '))
except ValueError:
    Fa = 10.0
try:
    af = float(input('anglar frequency for forced oscillation [1/s] (default=2.0): '))
except ValueError:
    af = 2.0

# initial condition
x0 = 0.0
v0 = 0.0
s0 = [x0, v0]               # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(forcedHarmonicOscillator, s0, t, args=(k,m,c,Fa,af))
xo = sol[:, 0]
xamp = np.max(xo)

xi = Fa*np.sin(af*t)     # 入力信号

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-5, 50), ylim=(-5, 5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
amp_i_template = r'$A_{{in}}$ = {0:.1f} N'.format(Fa)
amp_o_template = r'$A_{{out}}$ = {0:.1f} m'.format(xamp)
ax.text(0.75, 0.3, amp_i_template, transform=ax.transAxes)
ax.text(0.75, 0.5, amp_o_template, transform=ax.transAxes)

var_template = r'$k$ = {0:.1f} N/m, $m$ = {1:.1f} kg, $c$ = {2:.1f} kg/s'.format(k,m,c)
ax.text(0.4, 0.9, var_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

period_template = r'$\tau$ = {0:.2f} s, $\omega_r$ = {1:.2f} s$^{{-1}}$ ({2}), $\omega_f$ = {3:.2f} s$^{{-1}}$'.format(tau,afreq,cond,af)
ax.text(0.1, 0.8, period_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

input, = plt.plot([], [], 'bo-', animated=True)
output, = plt.plot([], [], 'ro-', animated=True)
# ここでは[],[]としているが、下でline.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

time_template = '$t$ = %.1f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return input, output, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_i = [0, l + xi[i]]
    y_i = [-2, -2]
    x_o = [0, l + xo[i]]
    y_o = [0, 0]
    input.set_data(x_i,y_i)
    output.set_data(x_o,y_o)
    time_text.set_text(time_template % (i*dt))
    return input, output, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/forced_harmonicOsci_(k={0:.1f},m={1:.1f},c={2:.1f},af={3:.2f}).gif'.format(k,m,c,af)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()