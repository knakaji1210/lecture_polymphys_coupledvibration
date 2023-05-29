# ordinary differential equation of damped ladder vibration (transverse)

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def dampedLadderVibration(s, t, k, m, c):
    [x1, v1, x2, v2, x3, v3, x4, v4, x5, v5, x6, v6, x7, v7, x8, v8, x9, v9]  = s
    a1 =        - 2*k*x1/m + k*x2/m - c*v1/m
    a2 = k*x1/m - 2*k*x2/m + k*x3/m - c*v2/m
    a3 = k*x2/m - 2*k*x3/m + k*x4/m - c*v3/m
    a4 = k*x3/m - 2*k*x4/m + k*x5/m - c*v4/m
    a5 = k*x4/m - 2*k*x5/m + k*x6/m - c*v5/m
    a6 = k*x5/m - 2*k*x6/m + k*x7/m - c*v6/m
    a7 = k*x6/m - 2*k*x7/m + k*x8/m - c*v7/m
    a8 = k*x7/m - 2*k*x8/m + k*x9/m - c*v8/m
    a9 = k*x8/m - 2*k*x9/m          - c*v9/m
    dsdt = [v1, a1, v2, a2, v3, a3, v4, a4, v5, a5, v6, a6, v7, a7, v8, a8, v9, a9]
    return dsdt

'''
上のladderVibrationをもっと簡潔に書けないかは未決
'''

# variables
try:
    k = float(input('spring constant [N/m] (default=100.0): '))
except ValueError:
    k = 100.0                # [N/m] spring constant
try:
    m = float(input('mass (default=20.0): '))    # [kg] mass
except ValueError:
    m = 20.0
try:
    c = float(input('damping coefficient (default=10.0): '))    # [kg/s] damping coefficient
except ValueError:
    c = 10.0
l = 10                      # [m] equilibrium length for each component
n = 9                       # number of mass
pos = [i*l for i in range(n+2)]
L = pos[-1]                 # [m] total length
af = [2*np.sqrt(k/m)*np.sin((i+1)*np.pi/(2*(n+1))) for i in range(n)]    # angular frequencies                            # periods of af [s] (T)
rho = c/(2*m)
tau = 1/rho
try:
    mode2show = int(input('mode to visualize (default=4): '))
except ValueError:
    mode2show = 4
if np.abs(rho - af[0]) < 0.05: # critical damping for q1 (rho == af[0]はほぼ無理)
    peri1 = tau                # period [s] (T)
    cond1 = "cd"
elif rho - af[0] < 0: # under damping for q1
    af1 = np.sqrt(af[0]**2 - rho**2)
    peri1 = 2*np.pi/af1      # period of af1[s] (T1)
    cond1 = "ud"
elif rho - af[0] > 0: # over damping for q1
    eta1 = np.sqrt(rho**2 - af[0]**2)
    peri1 = 1/(rho - eta1)   # period of af1[s] (T1)
    cond1 = "od"
if np.abs(rho - af[1]) < 0.05: # critical damping for q2 (rho == af[1]はほぼ無理)
    peri2 = tau                # period [s] (T)
    cond2 = "cd"
elif rho - af[1] < 0: # under damping for q1
    af2 = np.sqrt(af[1]**2 - rho**2)
    peri2 = 2*np.pi/af2      # period of af1[s] (T1)
    cond2 = "ud"
elif rho - af[1] > 0: # over damping for q1
    eta2 = np.sqrt(rho**2 - af[1]**2)
    peri2 = 1/(rho - eta2)   # period of af1[s] (T1)
    cond2 = "od"
if np.abs(rho - af[2]) < 0.05: # critical damping for q3 (rho == af[2]はほぼ無理)
    peri2 = tau                # period [s] (T)
    cond2 = "cd"
elif rho - af[2] < 0: # under damping for q1
    af3 = np.sqrt(af[2]**2 - rho**2)
    peri3 = 2*np.pi/af3      # period of af1[s] (T1)
    cond3 = "ud"
elif rho - af[2] > 0: # over damping for q1
    eta3 = np.sqrt(rho**2 - af[2]**2)
    peri3 = 1/(rho - eta3)   # period of af1[s] (T1)
    cond3 = "od"
if np.abs(rho - af[mode2show-1]) < 0.05: # critical damping for q_n (rho == af[n-1]はほぼ無理)
    perin = tau                # period [s] (T)
    condn = "cd"
elif rho - af[mode2show-1] < 0:
    afn = np.sqrt(af[mode2show-1]**2 - rho**2)
    perin = 2*np.pi/afn      # period of af_n[s] (T_n) 
    condn = "ud"
elif rho - af[mode2show-1] > 0:
    etan = np.sqrt(rho**2 - af[mode2show-1]**2)
    perin = 1/(rho - etan)   # period of af1[s] (T1) 
    condn = "od"

if cond1 == "cd" or cond2 == "cd" or cond3 == "cd" or condn == "cd":
    tmax = 8*peri1
else:
    tmax = 4*np.max([peri1,peri2,peri3,perin])  # [s] duration time   

dt = 0.05                   # [s] interval time

# initial condition
amp = 5
try:
    ini = str(input('initial condition: "s"ingle, "m"ultiple, "d"elta : '))
except ValueError:
    ini = "s"
if ini == "s":
    try:
        mode = int(input('mode to stimulate (default=1): '))
    except ValueError:
        mode = 1
    x_ini = [amp*np.sin(mode*np.pi*(i+1)/(n+1)) for i in range(n)]     # [m]   initial position
elif ini == "m":
    x_ini = [amp*np.sin(np.pi*(i+1)/(n+1)) + amp*np.sin(3*np.pi*(i+1)/(n+1)) + amp*np.sin(5*np.pi*(i+1)/(n+1)) for i in range(n)]     # [m]   initial position
elif ini == "d":
    x_ini = [2*amp]+np.zeros(n).tolist()     # [m]   initial position
else:
    x_ini = np.zeros(n+1).tolist()     # [m]   initial position
v_ini = np.zeros(n)                                             # [m/s]   initial velocity
xv_ini = [[x,v] for (x,v) in zip(x_ini,v_ini)]
s0 = [x for component in xv_ini for x in component]

t = np.arange(0, tmax, dt)

sol = odeint(dampedLadderVibration, s0, t, args=(k, m, c))
x1, x2, x3, x4, x5, x6, x7, x8, x9 = sol[:, 0], sol[:, 2], sol[:, 4], sol[:, 6], sol[:, 8], sol[:, 10], sol[:, 12], sol[:, 14], sol[:, 16]
q1 = np.sin(np.pi*1/10)*x1 + np.sin(np.pi*2/10)*x2 + np.sin(np.pi*3/10)*x3 + np.sin(np.pi*4/10)*x4 + np.sin(np.pi*5/10)*x5 + np.sin(np.pi*6/10)*x6 + np.sin(np.pi*7/10)*x7 + np.sin(np.pi*8/10)*x8 + np.sin(np.pi*9/10)*x9                      # normal mode q1
q2 = np.sin(2*np.pi*1/10)*x1 + np.sin(2*np.pi*2/10)*x2 + np.sin(2*np.pi*3/10)*x3 + np.sin(2*np.pi*4/10)*x4 + np.sin(2*np.pi*5/10)*x5 + np.sin(2*np.pi*6/10)*x6 + np.sin(2*np.pi*7/10)*x7 + np.sin(2*np.pi*8/10)*x8 + np.sin(2*np.pi*9/10)*x9    # normal mode q2
q3 = np.sin(3*np.pi*1/10)*x1 + np.sin(3*np.pi*2/10)*x2 + np.sin(3*np.pi*3/10)*x3 + np.sin(3*np.pi*4/10)*x4 + np.sin(3*np.pi*5/10)*x5 + np.sin(3*np.pi*6/10)*x6 + np.sin(3*np.pi*7/10)*x7 + np.sin(3*np.pi*8/10)*x8 + np.sin(3*np.pi*9/10)*x9    # normal mode q3
qn = np.sin(mode2show*np.pi*1/10)*x1 + np.sin(mode2show*np.pi*2/10)*x2 + np.sin(mode2show*np.pi*3/10)*x3 + np.sin(mode2show*np.pi*4/10)*x4 + np.sin(mode2show*np.pi*5/10)*x5 + np.sin(mode2show*np.pi*6/10)*x6 + np.sin(mode2show*np.pi*7/10)*x7 + np.sin(mode2show*np.pi*8/10)*x8 + np.sin(mode2show*np.pi*9/10)*x9    # normal mode qn

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(0, L), ylim=(-L/5, 0.7*L))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')

line, = plt.plot([], [], 'ro-', animated=True)
norm1, = plt.plot([], [], 'bo-', animated=True)
norm2, = plt.plot([], [], 'go-', animated=True)
norm3, = plt.plot([], [], 'yo-', animated=True)
norm, = plt.plot([], [], 'co-', animated=True)
# ここでは[],[]としているが、下でline.set_dataなどで実際の値を入れている

var_template = r'$k$ = {0:.1f} N/m, $m$ = {1:.1f} kg, $c$ = {2:.1f} kg/s'.format(k,m,c)
var_text = ax.text(0.4, 0.92, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

peri1_template = r'$\tau$ = {0:.2f} s, $T_1$ = {1:.2f} s ({2}), $T_2$ = {3:.2f} s ({4})'.format(tau,peri1,cond1,peri2,cond2)
peri1_text = ax.text(0.1, 0.85, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

peri2_template = r'$T_3$ = {0:.2f} s ({1}), $T_{2}$ = {3:.2f} s ({4})'.format(peri3,cond3,mode2show,perin,condn)
peri2_text = ax.text(0.1, 0.78, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.92, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    peri1_text.set_text('')
    peri2_text.set_text('')
    var_text.set_text('')
    return line, norm1, norm2, norm3, norm, time_text, peri1_text, peri2_text, var_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line.set_data(pos, 
            [0, x1[i], x2[i], x3[i], x4[i], x5[i], x6[i], x7[i], x8[i], x9[i], 0])
    norm1.set_data([0, L/2 + q1[i]], [10, 10])
    norm2.set_data([0, L/2 + q2[i]], [20, 20])
    norm3.set_data([0, L/2 + q3[i]], [30, 30])
    norm.set_data([0, L/2 + qn[i]], [40, 40])
    time_text.set_text(time_template % (i*dt))
    peri1_text.set_text(peri1_template)
    peri2_text.set_text(peri2_template)
    var_text.set_text(var_template)
    return line, norm1, norm2, norm3, norm, time_text, peri1_text, peri2_text, var_text

f = np.arange(0, len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)

if ini == "s": 
    savefile = './gif/damped_ladderdVibration_single_(ini{0},show{1},k={2:.1f},m={3:.1f}c={4:.1f}).gif'.format(mode,mode2show,k,m,c)
    ani.save(savefile, writer='pillow', fps=fps)
elif ini == "m":
    savefile = './gif/damped_ladderdVibration_multiple_(show{0},k={1:.1f},m={2:.1f}c={3:.1f}).gif'.format(mode2show,k,m,c)
    ani.save(savefile, writer='pillow', fps=fps)
elif ini == "d":
    savefile = './gif/damped_ladderdVibration_delta_(show{0},k={1:.1f},m={2:.1f}c={3:.1f}).gif'.format(mode2show,k,m,c)
    ani.save(savefile, writer='pillow', fps=fps)
else:
    pass

plt.show()