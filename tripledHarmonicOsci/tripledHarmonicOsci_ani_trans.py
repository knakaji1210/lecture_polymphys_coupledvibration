# ordinary differential equation of tripled harmonic oscillator (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def tripledharmonicOscillator(s, t, k, m):

    x1, v1, x2, v2, x3, v3 = s              # s = (x1, v1, x2, v2, x3, v3)
    a1 = -2*k*x1/m + k*x2/m
    a2 = k*x1/m - 2*k*x2/m + k*x3/m
    a3 = k*x2/m - 2*k*x3/m
    dsdt = [v1, a1, v2, a2, v3, a3]
    return dsdt

# variables
k = 30                      # [N/m] spring constant
m = 10                      # [kg] mass
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
l4 = 20                     # [m] equilibrium length
L = l1 + l2 + l3 + l4
n = 3                       # number of mass
af = [2*np.sqrt(k/m)*np.sin((i+1)*np.pi/(2*(n+1))) for i in range(n)]    # angular frequencies
peri = [2*np.pi/af[i] for i in range(n)]                                 # periods of af [s] (T)
tmax = 2*peri[0]             # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
try:
    x1_0 = float(input('initial position of mass1 (default=5.0): '))
except ValueError:
    x1_0 = 5.0
try:
    x2_0 = float(input('initial position of mass2 (default=10.0): '))
except ValueError:
    x2_0 = 10.0
try:
    x3_0 = float(input('initial position of mass3 (default=-5.0): '))
except ValueError:
    x3_0 = -5.0
try:
    v1_0 = float(input('initial velocity of mass1 (default=0.0): '))
except ValueError:
    v1_0 = 0.0
try:
    v2_0 = float(input('initial velocity of mass2 (default=0.0): '))
except ValueError:
    v2_0 = 0.0
try:
    v3_0 = float(input('initial velocity of mass3 (default=0.0): '))
except ValueError:
    v3_0 = 0.0

s0 = [x1_0, v1_0, x2_0, v2_0, x3_0, v3_0]   # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(tripledharmonicOscillator, s0, t, args=(k, m))  # ODEの解を求めている
x1, x2, x3 = sol[:, 0], sol[:, 2], sol[:, 4]    # [x1], [x2], [x3]が出てくる
q1 = (1/2)*(x1 + np.sqrt(2)*x2+x3)              # normal mode q1
q2 = (1/2)*(-np.sqrt(2)*x1 + np.sqrt(2)*x3)     # normal mode q2
q3 = (1/2)*(x1 - np.sqrt(2)*x2+x3)              # normal mode q3

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/4, L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
norm3, = plt.plot([], [], 'yo-', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

peri_template = '$T_1$ = {0:.2f} s, $T_2$ = {1:.2f} s, $T_3$ = {2:.2f} s'.format(peri[0],peri[1],peri[2])
peri_text = ax.text(0.1, 0.85, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.92, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    peri_text.set_text('')
    return line, norm1, norm2, norm3, time_text, peri_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1, l1 + l2, l1 + l2 + l3, L], [0, x1[i], x2[i], x3[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [l1, l1])
    norm2.set_data([0, L/2 + q2[i]], [2*l1, 2*l1])
    norm3.set_data([0, L/2 + q3[i]], [3*l1, 3*l1])
    time_text.set_text(time_template % (i*dt))
    peri_text.set_text(peri_template)
    return line, norm1, norm2, norm3, time_text, peri_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/tripleddHarmonicOsci_t_(x1={0:.1f},x2={1:.1f},x3={2:.1f}).gif'.format(x1_0,x2_0,x3_0)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()