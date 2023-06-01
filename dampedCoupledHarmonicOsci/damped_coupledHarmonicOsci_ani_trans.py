# ordinary differential equation of damped coupled harmonic oscillator (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedCoupledHarmonicOscillator(s, t, k1, k2, m, c):

    x1, v1, x2, v2 = s              # s = (x1, v1, x2, v2)
    a1 = - (k1+k2)*x1/m + k2*x2/m - c*v1/m  
    a2 = k2*x1/m - (k1+k2)*x2/m - c*v2/m
    dsdt = [v1, a1, v2, a2]
    return dsdt

# variables
try:
    k1 = float(input('spring constant 1 [N/m] (default=10.0): '))
except ValueError:
    k1 = 10.0               # [N/m] spring constant 1
try:
    k2 = float(input('spring constant 2 [N/m] (default=20.0): '))
except ValueError:
    k2 = 20.0               # [N/m] spring constant 2
try:
    m = float(input('mass [kg] (default=1.0): '))
except ValueError:
    m = 1.0                 # [kg] mass
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5                 # [m/s] damping coefficient
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
L = l1 + l2 + l3
af1_0 = np.sqrt(k1/m)            # natural　angular frequency
af2_0 = np.sqrt((k1+2*k2)/m)     # natural　angular frequency
rho = c/(2*m)
tau = 1/rho
if np.abs(rho - af1_0) < 0.05: # critical damping (rho == af1_0はほぼ無理)
    peri1 = tau                # period [s] (T)
    cond1 = "cd"
elif rho - af1_0 < 0: # under damping for q1
    af1 = np.sqrt(af1_0**2 - rho**2)
    peri1 = 2*np.pi/af1      # period of af1[s] (T1)
    cond1 = "ud"
elif rho - af1_0 > 0: # over damping for q1
    eta1 = np.sqrt(rho**2 - af1_0**2)
    peri1 = 1/(rho - eta1)      # period [s] (T)
    cond1 = "od"
if np.abs(rho - af2_0) < 0.05: # critical damping (rho == af2_0はほぼ無理)
    peri2 = tau                # period [s] (T)
    cond2 = "cd"
elif rho - af2_0 < 0: # under damping for q2
    af2 = np.sqrt(af2_0**2 - rho**2)
    peri2 = 2*np.pi/af2      # period of af2[s] (T2) 
    cond2 = "ud"
elif rho - af2_0 > 0:
    eta2 = np.sqrt(rho**2 - af2_0**2)
    peri2 = 1/(rho - eta2)      # period [s] (T)
    cond2 = "od"

if cond1 == "cd" or cond2 == "od":
    tmax = 8*peri1
else:
    tmax = 4*np.max([peri1,peri2])              # [s] duration time
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

sol = odeint(dampedCoupledHarmonicOscillator, s0, t, args=(k1, k2, m, c))  # ODEの解を求めている
x1, x2 = sol[:, 0], sol[:, 2]   # [x1], [x2]が出てくる
q1 = (1/np.sqrt(2))*(x1 + x2)   # normal mode q1
q2 = (1/np.sqrt(2))*(-x1 + x2)  # normal mode q2

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/4, 1.2*L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
var1_template = r'$k, \kappa$ = {0:.1f}, {1:.1f} N/m'.format(k1,k2)
ax.text(0.6, 0.92, var1_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
var2_template = r'$m$ = {0:.1f} kg, $c$ = {1:.1f} kg/s'.format(m,c)
ax.text(0.6, 0.85, var2_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
peri_template = r'$\tau$ = {0:.2f} s, $T_1$ = {1:.2f} s ({2}), $T_2$ = {3:.2f} s ({4})'.format(tau,peri1,cond1,peri2,cond2)
ax.text(0.1, 0.78, peri_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.92, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return line, norm1, norm2, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1, l1 + l2, L], [0, x1[i], x2[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [20, 20])
    norm2.set_data([0, L/2 + q2[i]], [40, 40])
    time_text.set_text(time_template % (i*dt))
    return line, norm1, norm2, time_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/damped_coupledHarmonicOsci_t_(x1={0:.1f},x2={1:.1f},k1={2:.1f},k2={3:.1f},m={4:.1f}c={5:.1f}).gif'.format(x1_0,x2_0,k1,k2,m,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()