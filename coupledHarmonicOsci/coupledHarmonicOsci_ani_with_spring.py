# ordinary differential equation of coupled harmonic oscillator

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
    k2 = float(input('spring constant 2 [N/m] (default=20.0): '))
except ValueError:
    k2 = 20.0               # [N/m] spring constant
try:
    m = float(input('mass [kg] (default=1.0): '))
except ValueError:
    m = 1.0                # [kg] mass
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
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-1, 3))
#ax = fig.add_subplot(111, xlim=(0, L), ylim=(-1.5, 1.5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
var1_template = r'$k, \kappa$ = {0:.1f}, {1:.1f} N/m'.format(k1,k2)
ax.text(0.6, 0.9, var1_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
var2_template = r'$m$ = {0:.1f} kg'.format(m)
ax.text(0.6, 0.8, var2_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
peri_template = '$T_1$ = {0:.2f} s, $T_2$ = {1:.2f} s'.format(peri1,peri2)
ax.text(0.1, 0.8, peri_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
x_rod_l, x_rod_r = [0, l1/4], [L-l3/4, L]   # 両端のロッド
y = [0, 0]
ax.plot(x_rod_l,y, c='b')
ax.plot(x_rod_r,y, c='b')

mass, = plt.plot([], [], 'ro', markersize='10', animated=True, zorder=10)
rod1, = ax.plot([],[], 'b', animated=True)
rod2, = ax.plot([],[], 'b', animated=True)
tri1, = ax.plot([],[], 'b', animated=True)
tri2, = ax.plot([],[], 'b', animated=True)
tri3, = ax.plot([],[], 'b', animated=True)
norm1, = plt.plot([], [], 'bo-', markersize='10', animated=True)
norm2, = plt.plot([], [], 'go-', markersize='10', animated=True)
# ここでは[],[]としているが、下でlinei.set_dataで実際の値を入れている

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

# 基準モード描画あり
def init_w():               # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return mass, rod1, rod2, tri1, tri2, tri3, norm1, norm2, time_text

# 基準モード描画なし
def init_wo():              # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return mass, rod1, rod2, tri1, tri2, tri3, time_text

# 基準モード描画あり
def update_w(i):             # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod1 = [3*l1/4 + x1[i], l1 + x1[i]+l2/4]
    x_rod2 = [l1 + l2 + x2[i]-l2/4, l1 + l2 + x2[i]+l3/4]
    rod1.set_data(x_rod1,y)
    rod2.set_data(x_rod2,y)
    x_tri1 = np.linspace(l1/4, 3*l1/4 + x1[i],100)
    y_tri1 = (1/6)*np.arccos(np.cos(6*np.pi*(x_tri1 - l1/4)/(x1[i]+l1/2)-np.pi/2+0.1))-1/4
    tri1.set_data(x_tri1,y_tri1)
    x_tri2 = np.linspace(l1 + x1[i]+l2/4, l1 + 3*l2/4 + x2[i],100)
    y_tri2 = (1/6)*np.arccos(np.cos(6*np.pi*(x_tri2 - l1 - x1[i] - l2/4)/(x2[i]-x1[i]+l2/2)-np.pi/2+0.1))-1/4
    tri2.set_data(x_tri2,y_tri2)
    x_tri3 = np.linspace(l1 + l2 + x2[i]+l3/4, L-l3/4,100)
    y_tri3 = (1/6)*np.arccos(np.cos(6*np.pi*(x_tri3 - l1 - l2 - x2[i] - l3/4)/(l3/2 - x2[i])-np.pi/2+0.1))-1/4
    tri3.set_data(x_tri3,y_tri3)
    mass.set_data([0, l1 + x1[i], l1 + l2 + x2[i], L], [0, 0, 0, 0])
    norm1.set_data([0, L/2 + q1[i]], [1, 1])    # L/2を中心に描画
    norm2.set_data([0, L/2 + q2[i]], [2, 2])
    time_text.set_text(time_template % (i*dt))
    return mass, rod1, rod2, tri1, tri2, tri3, norm1, norm2, time_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

# 基準モード描画なし
def update_wo(i):            # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    x_rod1 = [3*l1/4 + x1[i], l1 + x1[i]+l2/4]
    x_rod2 = [l1 + l2 + x2[i]-l2/4, l1 + l2 + x2[i]+l3/4]
    rod1.set_data(x_rod1,y)
    rod2.set_data(x_rod2,y)
    x_tri1 = np.linspace(l1/4, 3*l1/4 + x1[i],100)
    y_tri1 = (1/6)*np.arccos(np.cos(6*np.pi*(x_tri1 - l1/4)/(x1[i]+l1/2)-np.pi/2+0.1))-1/4
    tri1.set_data(x_tri1,y_tri1)
    x_tri2 = np.linspace(l1 + x1[i]+l2/4, l1 + 3*l2/4 + x2[i],100)
    y_tri2 = (1/6)*np.arccos(np.cos(6*np.pi*(x_tri2 - l1 - x1[i] - l2/4)/(x2[i]-x1[i]+l2/2)-np.pi/2+0.1))-1/4
    tri2.set_data(x_tri2,y_tri2)
    x_tri3 = np.linspace(l1 + l2 + x2[i]+l3/4, L-l3/4,100)
    y_tri3 = (1/6)*np.arccos(np.cos(6*np.pi*(x_tri3 - l1 - l2 - x2[i] - l3/4)/(l3/2 - x2[i])-np.pi/2+0.1))-1/4
    tri3.set_data(x_tri3,y_tri3)
    mass.set_data([0, l1 + x1[i], l1 + l2 + x2[i], L], [0, 0, 0, 0])
    time_text.set_text(time_template % (i*dt))
    return mass, rod1, rod2, tri1, tri2, tri3, time_text

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
    savefile = './gif/coupledHarmonicOsci_wn_wsp_(x1={0:.1f},v1={1:.1f},x2={2:.1f},v2={3:.1f},k1={4:.1f},k2={5:.1f},m={6:.1f}).gif'.format(x1_0,v1_0,x2_0,v2_0,k1,k2,m)
    ani.save(savefile, writer='pillow', fps=fps)
    plt.show()

elif n_mode == "n": # 基準モードを描画しない場合
    ani = FuncAnimation(fig, update_wo, frames=f,
                    init_func=init_wo, blit=True, interval=frame_int, repeat=True)
    savefile = './gif/coupledHarmonicOsci_won_wsp_(x1={0:.1f},v1={1:.1f},x2={2:.1f},v2={3:.1f},k1={4:.1f},k2={5:.1f},m={6:.1f}).gif'.format(x1_0,v1_0,x2_0,v2_0,k1,k2,m)
    ani.save(savefile, writer='pillow', fps=fps)
    plt.show()
else:
    print('Select "y" or "n"!')
    pass