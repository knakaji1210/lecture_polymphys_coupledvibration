# ordinary differential equation of damped coupled harmonic oscillator

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
k1 = 100                    # [N/m] spring constant
k2 = 200                    # [N/m] spring constant
m = 20                      # [kg] mass
try:
    c = float(input('damping coefficient (default=10.0): '))    # [kg/s] damping coefficient
except ValueError:
    c = 10.0
l1 = 20                     # [m] equilibrium length
l2 = 20                     # [m] equilibrium length
l3 = 20                     # [m] equilibrium length
L = l1 + l2 + l3
af1_0 = np.sqrt(k1/m)            # natural　angular frequency
af2_0 = np.sqrt((k1+2*k2)/m)     # natural　angular frequency
rho = c/(2*m)
tau = 1/rho
if rho - af1_0 < 0:
    af1 = np.sqrt(af1_0**2 - rho**2)
    peri1 = 2*np.pi/af1      # period of af1[s] (T1)
    cond1 = "underdamped"
elif rho - af1_0 > 0:
    af1 = np.sqrt(rho**2 - af1_0**2)
    peri1 = 2*np.pi/af1      # period of af1[s] (T1)
    cond1 = "overdamped"
if rho - af2_0 < 0:
    af2 = np.sqrt(af2_0**2 - rho**2)
    peri2 = 2*np.pi/af2      # period of af2[s] (T2) 
    cond2 = "underdamped"
elif rho - af2_0 > 0:
    af2 = np.sqrt(rho**2 - af2_0**2)
    peri2 = 2*np.pi/af2      # period of af2[s] (T2) 
    cond2 = "overdamped"
'''
とりあえずrho = afreq0（臨界減衰）のときは無視
'''
tmax = 4*peri1              # [s] duration time
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
    v1_0 = float(input('initial velocity of mass1 (default=0.0): '))
except ValueError:
    v1_0 = 0.0
try:
    v2_0 = float(input('initial velocity of mass2 (default=0.0): '))
except ValueError:
    v2_0 = 0.0

s0 = [x1_0, v1_0, x2_0, v2_0]   # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(dampedCoupledHarmonicOscillator, s0, t, args=(k1, k2, m, c))  # ODEの解を求めている
x1, x2 = sol[:, 0], sol[:, 2]   # [x1], [x2]が出てくる
q1 = (1/np.sqrt(2))*(x1 + x2)   # normal mode q1
q2 = (1/np.sqrt(2))*(-x1 + x2)  # normal mode q2

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-1, 3.5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

damp_template = r'$c$ = {0:.1f} kg/s, $\tau$ = {1:.2f} s'.format(c,tau)
damp_text = ax.text(0.1, 0.8, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

peri_template = r'$T_1$ = {0:.2f} s ({1}), $T_2$ = {2:.2f} s ({3})'.format(peri1,cond1,peri2,cond2)
peri_text = ax.text(0.1, 0.7, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# 基準モード描画あり
def init_w():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    damp_text.set_text('')
    peri_text.set_text('')
    return line, norm1, norm2, time_text, damp_text, peri_text

# 基準モード描画なし
def init_wo():              # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    damp_text.set_text('')
    peri_text.set_text('')
    return line, time_text, damp_text, peri_text

# 基準モード描画あり
def update_w(i):             # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1 + x1[i], l1 + l2 + x2[i], L], [0, 0, 0, 0])
    norm1.set_data([0, L/2 + q1[i]], [1, 1])    # L/2を中心に描画
    norm2.set_data([0, L/2 + q2[i]], [2, 2])
    time_text.set_text(time_template % (i*dt))
    damp_text.set_text(damp_template)
    peri_text.set_text(peri_template)
    return line, norm1, norm2, time_text, damp_text, peri_text

# 基準モード描画なし
def update_wo(i):            # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, l1 + x1[i], l1 + l2 + x2[i], L], [0, 0, 0, 0])
    time_text.set_text(time_template % (i*dt))
    damp_text.set_text(damp_template)
    peri_text.set_text(peri_template)
    return line, time_text, damp_text, peri_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

try:
    n_mode = input('With normal modes? (y or n): ')
except:
    n_mode = "y"

if n_mode == "y": # 基準モードを描画する場合
    ani = FuncAnimation(fig, update_w, frames=f,
                    init_func=init_w, blit=True, interval=frame_int, repeat=True)
    savefile = './gif/damped_coupledHarmonicOsci_wn_(x1={0:.1f},v1={1:.1f},x2={2:.1f},v2={3:.1f},c={4:.1f}).gif'.format(x1_0,v1_0,x2_0,v2_0,c)
    ani.save(savefile, writer='pillow', fps=fps)
    plt.show()

elif n_mode == "n": # 基準モードを描画しない場合
    ani = FuncAnimation(fig, update_wo, frames=f,
                    init_func=init_wo, blit=True, interval=frame_int, repeat=True)
    savefile = './gif/damped_coupledHarmonicOsci_won_(x1={0:.1f},v1={1:.1f},x2={2:.1f},v2={3:.1f},c={4:.1f}).gif'.format(x1_0,v1_0,x2_0,v2_0,c)
    ani.save(savefile, writer='pillow', fps=fps)
    plt.show()
else:
    print('Select "y" or "n"!')
    pass

