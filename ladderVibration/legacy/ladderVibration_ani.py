# ordinary differential equation of ladder vibration (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def ladderVibration(s, t, k, m):

    x1, v1, x2, v2, x3, v3, x4, v4, x5, v5, x6, v6, x7, v7, x8, v8, x9, v9  = s
    a1 =        - 2*k*x1/m + k*x2/m
    a2 = k*x1/m - 2*k*x2/m + k*x3/m
    a3 = k*x2/m - 2*k*x3/m + k*x4/m
    a4 = k*x3/m - 2*k*x4/m + k*x5/m
    a5 = k*x4/m - 2*k*x5/m + k*x6/m
    a6 = k*x5/m - 2*k*x6/m + k*x7/m
    a7 = k*x6/m - 2*k*x7/m + k*x8/m
    a8 = k*x7/m - 2*k*x8/m + k*x9/m
    a9 = k*x8/m - 2*k*x9/m        
    dsdt = [v1, a1, v2, a2, v3, a3, v4, a4, v5, a5, v6, a6, v7, a7, v8, a8, v9, a9]
    return dsdt

# variables
k = 60                      # [N/m] spring constant
m = 10                      # [kg] mass
l1 = 10                     # [m] equilibrium length
l2 = 10                     # [m] equilibrium length
l3 = 10                     # [m] equilibrium length
l4 = 10                     # [m] equilibrium length
l5 = 10                     # [m] equilibrium length
l6 = 10                     # [m] equilibrium length
l7 = 10                     # [m] equilibrium length
l8 = 10                     # [m] equilibrium length
l9 = 10                     # [m] equilibrium length
l10 = 10                    # [m] equilibrium length
L = l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9 + l10
n = 9                       # number of mass
af = [2*np.sqrt(k/m)*np.sin((i+1)*np.pi/(2*(n+1))) for i in range(n)]    # angular frequencies
peri = [2*np.pi/af[i] for i in range(n)]                                 # periods of af [s] (T)
tmax = 2*peri[0]            # [s] duration time
dt = 0.05                   # [s] interval time

# initial condition
amp = 5
mode = 2
x1_0 = amp*np.sin(mode*np.pi*1/10)   # [m]   initial position
x2_0 = amp*np.sin(mode*np.pi*2/10)   # [m]   initial position
x3_0 = amp*np.sin(mode*np.pi*3/10)   # [m]   initial position
x4_0 = amp*np.sin(mode*np.pi*4/10)   # [m]   initial position
x5_0 = amp*np.sin(mode*np.pi*5/10)   # [m]   initial position
x6_0 = amp*np.sin(mode*np.pi*6/10)   # [m]   initial position
x7_0 = amp*np.sin(mode*np.pi*7/10)   # [m]   initial position
x8_0 = amp*np.sin(mode*np.pi*8/10)   # [m]   initial position
x9_0 = amp*np.sin(mode*np.pi*9/10)   # [m]   initial position

#x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0 = 10, 0, 0, 0, 0, 0, 0, 0, 0

#x1_0 = amp*np.sin(mode*np.pi*1/10) + amp*np.sin(3*np.pi*1/10) + amp*np.sin(5*np.pi*1/10)    # [m]   initial position
#x2_0 = amp*np.sin(mode*np.pi*2/10) + amp*np.sin(3*np.pi*2/10) + amp*np.sin(5*np.pi*2/10)     # [m]   initial position
#x3_0 = amp*np.sin(mode*np.pi*3/10) + amp*np.sin(3*np.pi*3/10) + amp*np.sin(5*np.pi*3/10)     # [m]   initial position
#x4_0 = amp*np.sin(mode*np.pi*4/10) + amp*np.sin(3*np.pi*4/10) + amp*np.sin(5*np.pi*4/10)     # [m]   initial position
#x5_0 = amp*np.sin(mode*np.pi*5/10) + amp*np.sin(3*np.pi*5/10) + amp*np.sin(5*np.pi*5/10)     # [m]   initial position
#x6_0 = amp*np.sin(mode*np.pi*6/10) + amp*np.sin(3*np.pi*6/10) + amp*np.sin(5*np.pi*6/10)     # [m]   initial position
#x7_0 = amp*np.sin(mode*np.pi*7/10) + amp*np.sin(3*np.pi*7/10) + amp*np.sin(5*np.pi*7/10)     # [m]   initial position
#x8_0 = amp*np.sin(mode*np.pi*8/10) + amp*np.sin(3*np.pi*8/10) + amp*np.sin(5*np.pi*8/10)     # [m]   initial position
#x9_0 = amp*np.sin(mode*np.pi*9/10) + amp*np.sin(3*np.pi*9/10) + amp*np.sin(5*np.pi*9/10)     # [m]   initial position

v1_0 = 0                    # [m/s] initial velocity
v2_0 = 0                    # [m/s] initial velocity
v3_0 = 0                    # [m/s] initial velocity
v4_0 = 0                    # [m/s] initial velocity
v5_0 = 0                    # [m/s] initial velocity
v6_0 = 0                    # [m/s] initial velocity
v7_0 = 0                    # [m/s] initial velocity
v8_0 = 0                    # [m/s] initial velocity
v9_0 = 0                    # [m/s] initial velocity

s0 = [x1_0, v1_0, x2_0, v2_0, x3_0, v3_0, x4_0, v4_0, x5_0, v5_0, x6_0, v6_0, x7_0, v7_0, x8_0, v8_0, x9_0, v9_0]   # initial condition

t = np.arange(0, tmax, dt)

sol = odeint(ladderVibration, s0, t, args=(k, m))
x1, x2, x3, x4, x5, x6, x7, x8, x9 = sol[:, 0], sol[:, 2], sol[:, 4], sol[:, 6], sol[:, 8], sol[:, 10], sol[:, 12], sol[:, 14], sol[:, 16]
q1 = np.sin(np.pi*1/10)*x1 + np.sin(np.pi*2/10)*x2 + np.sin(np.pi*3/10)*x3 + np.sin(np.pi*4/10)*x4 + np.sin(np.pi*5/10)*x5 + np.sin(np.pi*6/10)*x6 + np.sin(np.pi*7/10)*x7 + np.sin(np.pi*8/10)*x8 + np.sin(np.pi*9/10)*x9                      # normal mode q1
q2 = np.sin(2*np.pi*1/10)*x1 + np.sin(2*np.pi*2/10)*x2 + np.sin(2*np.pi*3/10)*x3 + np.sin(2*np.pi*4/10)*x4 + np.sin(2*np.pi*5/10)*x5 + np.sin(2*np.pi*6/10)*x6 + np.sin(2*np.pi*7/10)*x7 + np.sin(2*np.pi*8/10)*x8 + np.sin(2*np.pi*9/10)*x9    # normal mode q2
q3 = np.sin(3*np.pi*1/10)*x1 + np.sin(3*np.pi*2/10)*x2 + np.sin(3*np.pi*3/10)*x3 + np.sin(3*np.pi*4/10)*x4 + np.sin(3*np.pi*5/10)*x5 + np.sin(3*np.pi*6/10)*x6 + np.sin(3*np.pi*7/10)*x7 + np.sin(3*np.pi*8/10)*x8 + np.sin(3*np.pi*9/10)*x9    # normal mode q3
q4 = np.sin(4*np.pi*1/10)*x1 + np.sin(4*np.pi*2/10)*x2 + np.sin(4*np.pi*3/10)*x3 + np.sin(4*np.pi*4/10)*x4 + np.sin(4*np.pi*5/10)*x5 + np.sin(4*np.pi*6/10)*x6 + np.sin(4*np.pi*7/10)*x7 + np.sin(4*np.pi*8/10)*x8 + np.sin(4*np.pi*9/10)*x9    # normal mode q_n

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/5, L/2))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
norm3, = plt.plot([], [], 'yo-', animated=True)
norm, = plt.plot([], [], 'co-', animated=True)
# ここでは[],[]としているが、下でline.set_dataなどで実際の値を入れている

peri_template = '$T_1$ = {0:.2f} s, $T_2$ = {1:.2f} s, $T_3$ = {2:.2f} s, $T_4$ = {3:.2f} s'.format(peri[0],peri[1],peri[2],peri[3])
peri_text = ax.text(0.1, 0.88, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.95, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    peri_text.set_text('')
    return line, norm1, norm2, norm3, norm, time_text, peri_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data([0, 
            l1, 
            l1 + l2, 
            l1 + l2 + l3, 
            l1 + l2 + l3 +l4, 
            l1 + l2 + l3 +l4 + l5,
            l1 + l2 + l3 +l4 + l5 + l6, 
            l1 + l2 + l3 +l4 + l5 + l6 + l7, 
            l1 + l2 + l3 +l4 + l5 + l6 + l7 + l8,
            l1 + l2 + l3 +l4 + l5 + l6 + l7 + l8 + l9,
             L], 
            [0, x1[i], x2[i], x3[i], x4[i], x5[i], x6[i], x7[i], x8[i], x9[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [10, 10])
    norm2.set_data([0, L/2 + q2[i]], [20, 20])
    norm3.set_data([0, L/2 + q3[i]], [30, 30])
    norm.set_data([0, L/2 + q4[i]], [40, 40])
    time_text.set_text(time_template % (i*dt))
    peri_text.set_text(peri_template)
    return line, norm1, norm2, norm3, norm, time_text, peri_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/ladderdVibration.gif'
ani.save(savefile, writer='pillow', fps=fps)

plt.show()