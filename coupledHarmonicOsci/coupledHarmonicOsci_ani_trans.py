# ordinary differential equation of coupled harmonic oscillator (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def coupledHarmonicOscillator(s, t, k1, k2, m):

    x1, v1, x2, v2 = s              # s = (x1, v1, x2, v2)
    a1 = - (k1+k2)*x1/m + k2*x2/m   
    a2 = k2*x1/m - (k1+k2)*x2/m
    dsdt = [v1, a1, v2, a2]
    return dsdt

# variables
try:
    k1 = float(input('spring constant 1 [N/m] (default=10.0): '))
except ValueError:
    k1 = 10.0               # [N/m] spring constant
try:
    k2 = float(input('spring constant 2 [N/m] (default=50.0): '))
except ValueError:
    k2 = 50.0               # [N/m] spring constant
try:
    m = float(input('mass [kg] (default=10.0): '))
except ValueError:
    m = 10.0                # [kg] mass
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
L = l1 + l2 + l3
af1 = np.sqrt(k1/m)         # angular frequency   
af2 = np.sqrt((k1+2*k2)/m)  # angular frequency
peri1 = 2*np.pi/af1         # period of af1[s] (T1)
peri2 = 2*np.pi/af2         # period of af1[s] (T2)     
tmax = 2*peri1              # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
try:
    x1_0 = float(input('initial position of mass1 [m] (default=5.0): '))
except ValueError:
    x1_0 = 5.0
try:
    x2_0 = float(input('initial position of mass2 [m] (default=10.0): '))
except ValueError:
    x2_0 = 10.0
try:
    v1_0 = float(input('initial velocity of mass1 [m/s] (default=0.0): '))
except ValueError:
    v1_0 = 0.0
try:
    v2_0 = float(input('initial velocity of mass2 [m/s] (default=0.0): '))
except ValueError:
    v2_0 = 0.0

s0 = [x1_0, v1_0, x2_0, v2_0]   # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(coupledHarmonicOscillator, s0, t, args=(k1, k2, m))  # ODEの解を求めている
x1, x2 = sol[:, 0], sol[:, 2]   # [x1], [x2]が出てくる
q1 = (1/np.sqrt(2))*(x1 + x2)   # normal mode q1
q2 = (1/np.sqrt(2))*(-x1 + x2)  # normal mode q2

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/4, L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

var1_template = r'$k_1, k_2$ = {0:.1f}, {1:.1f} N/m'.format(k1,k2)
var1_text = ax.text(0.6, 0.9, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

var2_template = r'$m$ = {0:.1f} kg'.format(m)
var2_text = ax.text(0.6, 0.8, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

peri_template = '$T_1$ = {0:.2f} s, $T_2$ = {1:.2f} s'.format(peri1,peri2)
peri_text = ax.text(0.1, 0.8, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    peri_text.set_text('')
    var1_text.set_text('')
    var2_text.set_text('')
    return line, norm1, norm2, time_text, peri_text, var1_text, var2_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1, l1 + l2, L], [0, x1[i], x2[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [20, 20])
    norm2.set_data([0, L/2 + q2[i]], [40, 40])
    time_text.set_text(time_template % (i*dt))
    peri_text.set_text(peri_template)
    var1_text.set_text(var1_template)
    var2_text.set_text(var2_template)
    return line, norm1, norm2, time_text, peri_text, var1_text, var2_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/coupledHarmonicOsci_t_(x1={0:.1f},x2={1:.1f},k1={2:.1f},k2={3:.1f},m={4:.1f}).gif'.format(x1_0,x2_0,k1,k2,m)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()