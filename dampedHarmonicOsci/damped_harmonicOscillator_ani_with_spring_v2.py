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
k = 100                     # [N/m] spring constant
m = 20                      # [kg] mass
try:
    c = float(input('damping coefficient (default=10.0): '))    # damping coefficient
except ValueError:
    c = 10.0
l = 20                      # [m] equilibrium length
afreq0 = np.sqrt(k/m)       # natural　angular frequency
rho = c/(2*m)
tau = 1/rho
if rho - afreq0 < 0:
    afreq = np.sqrt(afreq0**2 - rho**2)
    period = 2*np.pi/afreq      # period [s] (T)
    cond = "underdamped"
elif rho - afreq0 > 0:
    afreq = np.sqrt(rho**2 - afreq0**2)
    period = 2*np.pi/afreq      # period [s] (T)
    cond = "overdamped"
'''
とりあえずrho = afreq0（臨界減衰）のときは無視
'''

tmax = 4*period             # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
try:
    x0 = float(input('initial position (default=5.0): '))
except ValueError:
    x0 = 5.0
try:
    v0 = float(input('initial velocity (default=0.0): '))
except ValueError:
    v0 = 0.0

ic = [x0, v0]               # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(dampedHarmonicOscillator, ic, t, args=(k,m,c))  # ODEの解を求めている
#print(sol.shape) # (len(t), 2)が出てくる。sol = [[x],[v]]ということ
x = sol[:, 0]    # [x]が出てくる
#print(x.shape)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-5, 50), ylim=(-5, 5))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
x_rod1 = [0, l/4]
y = [0, 0]
ax.plot(x_rod1,y, c='b')

rod, = ax.plot([],[], 'b', animated=True)
triangle, = ax.plot([],[], 'b', animated=True)
mass, = plt.plot([], [], 'ro', markersize='10', animated=True)
# ここでは[],[]としているが、下で***.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

period_template = r'$\tau$ = {0:.2f} s ({2}), $T$ = {1:.2f} s'.format(tau,period,cond)
period_text = ax.text(0.1, 0.8, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    period_text.set_text('')
    return rod, triangle, mass, time_text, period_text

def update(i):              # ここのiは下のframes=fに対応した引数になっている
    x_rod2 = [3*l/4 + x[i], l + x[i]]
    x_mass = [0, l + x[i]]
    rod.set_data(x_rod2,y)
    x_tri = np.linspace(l/4, 3*l/4 + x[i],100)
    y_tri = (2/3)*np.arccos(np.cos(6*np.pi*(x_tri - l/4)/(x[i]+l/2)-np.pi/2+0.1))-1
    triangle.set_data(x_tri,y_tri)
    mass.set_data(x_mass,y)
    time_text.set_text(time_template % (i*dt))
    period_text.set_text(period_template)
    return rod, triangle, mass, time_text, period_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/dampedharmonicOsci_wsp_(x0={0:.1f},v0={1:.1f},c={2:.1f}).gif'.format(x0,v0,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()