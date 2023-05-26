# ordinary differential equation of damped tripled harmonic oscillator (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedTripledHarmonicOscillator(s, t, k, m, c):

    x1, v1, x2, v2, x3, v3 = s              # s = (x1, v1, x2, v2, x3, v3)
    a1 = -2*k*x1/m + k*x2/m - c*v1/m
    a2 = k*x1/m - 2*k*x2/m + k*x3/m - c*v2/m
    a3 = k*x2/m - 2*k*x3/m - c*v3/m
    dsdt = [v1, a1, v2, a2, v3, a3]
    return dsdt

# variables
k = 100                     # [N/m] spring constant
m = 20                      # [kg] mass
try:
    c = float(input('damping coefficient (default=10.0): '))    # [kg/s] damping coefficient
except ValueError:
    c = 10.0
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
l4 = 20                     # [m] equilibrium length
L = l1 + l2 + l3 + l4
n = 3                       # number of mass
af = [2*np.sqrt(k/m)*np.sin((i+1)*np.pi/(2*(n+1))) for i in range(n)]    # angular frequencies
peri = [2*np.pi/af[i] for i in range(n)]                                 # periods of af [s] (T)
tmax = 4*peri[0]             # [s] duration time
rho = c/(2*m)
tau = 1/rho
if rho - af[0] < 0:
    af1 = np.sqrt(af[0]**2 - rho**2)
    peri1 = 2*np.pi/af1      # period of af1[s] (T1)
    cond1 = "under"
elif rho - af[0] > 0:
    af1 = np.sqrt(rho**2 - af[0]**2)
    peri1 = 2*np.pi/af1      # period of af1[s] (T1)
    cond1 = "over"
if rho - af[1] < 0:
    af2 = np.sqrt(af[1]**2 - rho**2)
    peri2 = 2*np.pi/af2      # period of af2[s] (T2) 
    cond2 = "under"
elif rho - af[1] > 0:
    af2 = np.sqrt(rho**2 - af[1]**2)
    peri2 = 2*np.pi/af2      # period of af2[s] (T2) 
    cond2 = "over"
if rho - af[2] < 0:
    af3 = np.sqrt(af[2]**2 - rho**2)
    peri3 = 2*np.pi/af3      # period of af2[s] (T2) 
    cond3 = "under"
elif rho - af[2] > 0:
    af3 = np.sqrt(rho**2 - af[2]**2)
    peri3 = 2*np.pi/af3      # period of af2[s] (T2) 
    cond3 = "over"
'''
とりあえずrho = afreq0（臨界減衰）のときは無視
'''
tmax = 4*np.max([peri1,peri2,peri3])    # [s] duration time
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

sol = odeint(dampedTripledHarmonicOscillator, s0, t, args=(k, m, c))  # ODEの解を求めている
x1, x2, x3 = sol[:, 0], sol[:, 2], sol[:, 4]    # [x1], [x2], [x3]が出てくる
q1 = (1/2)*(x1 + np.sqrt(2)*x2+x3)              # normal mode q1
q2 = (1/2)*(-np.sqrt(2)*x1 + np.sqrt(2)*x3)     # normal mode q2
q3 = (1/2)*(x1 - np.sqrt(2)*x2+x3)              # normal mode q3

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/4, 1.4*L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
norm3, = plt.plot([], [], 'yo-', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

damp_template = r'$c$ = {0:.1f} kg/s, $\tau$ = {1:.2f} s'.format(c,tau)
damp_text = ax.text(0.1, 0.8, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

peri_template = r'$T_1$ = {0:.2f} s ({1}), $T_2$ = {2:.2f} s ({3}), $T_3$ = {4:.2f} s ({5})'.format(peri1,cond1,peri2,cond2,peri3,cond3)
peri_text = ax.text(0.1, 0.7, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    damp_text.set_text('')
    peri_text.set_text('')
    return line, norm1, norm2, norm3, time_text, damp_text, peri_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1, l1 + l2, l1 + l2 + l3, L], [0, x1[i], x2[i], x3[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [l1, l1])
    norm2.set_data([0, L/2 + q2[i]], [2*l1, 2*l1])
    norm3.set_data([0, L/2 + q3[i]], [3*l1, 3*l1])
    time_text.set_text(time_template % (i*dt))
    damp_text.set_text(damp_template)
    peri_text.set_text(peri_template)
    return line, norm1, norm2, norm3, time_text, damp_text, peri_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/damped_tripledHarmonicOsci_t_(x1={0:.1f},x2={1:.1f},x3={2:.1f},c={3:.1f}).gif'.format(x1_0,x2_0,x3_0,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()