import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k, m):

    x1, v1, x2, v2, x3, v3, x4, v4, x5, v5, x6, v6, x7, v7, x8, v8, x9, v9  = s
    a1 = -2*k*x1/m + k*x2/m         - c*v1/m
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

k = 30                      # [N/m] spring constant
m = 0.01                      # [kg] mass
c = 5                       # air resistance
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
amp = 5
mode = 1
#x1_0 = amp*np.sin(mode*np.pi*1/10) + amp*np.sin(3*np.pi*1/10) + amp*np.sin(5*np.pi*1/10)    # [m]   initial position
#x2_0 = amp*np.sin(mode*np.pi*2/10) + amp*np.sin(3*np.pi*2/10) + amp*np.sin(5*np.pi*2/10)     # [m]   initial position
#x3_0 = amp*np.sin(mode*np.pi*3/10) + amp*np.sin(3*np.pi*3/10) + amp*np.sin(5*np.pi*3/10)     # [m]   initial position
#x4_0 = amp*np.sin(mode*np.pi*4/10) + amp*np.sin(3*np.pi*4/10) + amp*np.sin(5*np.pi*4/10)     # [m]   initial position
#x5_0 = amp*np.sin(mode*np.pi*5/10) + amp*np.sin(3*np.pi*5/10) + amp*np.sin(5*np.pi*5/10)     # [m]   initial position
#x6_0 = amp*np.sin(mode*np.pi*6/10) + amp*np.sin(3*np.pi*6/10) + amp*np.sin(5*np.pi*6/10)     # [m]   initial position
#x7_0 = amp*np.sin(mode*np.pi*7/10) + amp*np.sin(3*np.pi*7/10) + amp*np.sin(5*np.pi*7/10)     # [m]   initial position
#x8_0 = amp*np.sin(mode*np.pi*8/10) + amp*np.sin(3*np.pi*8/10) + amp*np.sin(5*np.pi*8/10)     # [m]   initial position
#x9_0 = amp*np.sin(mode*np.pi*9/10) + amp*np.sin(3*np.pi*9/10) + amp*np.sin(5*np.pi*9/10)     # [m]   initial position

x1_0 = amp*np.sin(mode*np.pi*1/10)
x2_0 = amp*np.sin(mode*np.pi*2/10)
x3_0 = amp*np.sin(mode*np.pi*3/10)
x4_0 = amp*np.sin(mode*np.pi*4/10)
x5_0 = amp*np.sin(mode*np.pi*5/10)
x6_0 = amp*np.sin(mode*np.pi*6/10)
x7_0 = amp*np.sin(mode*np.pi*7/10)
x8_0 = amp*np.sin(mode*np.pi*8/10)
x9_0 = amp*np.sin(mode*np.pi*9/10)

#x1_0, x2_0, x3_0, x4_0, x5_0, x6_0, x7_0, x8_0, x9_0 = 10, 0, 0, 0, 0, 0, 0, 0, 0

v1_0 = 0                    # [m/s] initial velocity
v2_0 = 0                    # [m/s] initial velocity
v3_0 = 0                    # [m/s] initial velocity
v4_0 = 0                    # [m/s] initial velocity
v5_0 = 0                    # [m/s] initial velocity
v6_0 = 0                    # [m/s] initial velocity
v7_0 = 0                    # [m/s] initial velocity
v8_0 = 0                    # [m/s] initial velocity
v9_0 = 0                    # [m/s] initial velocity
s0 = [x1_0, v1_0, x2_0, v2_0, x3_0, v3_0, x4_0, v4_0, x5_0, v5_0, x6_0, v6_0, x7_0, v7_0, x8_0, v8_0, x9_0, v9_0]
af = 2*np.sqrt(k/m)*np.sin(np.pi/20)    
t = 16*np.pi/af             # [s] duration time
dt = 0.01                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(harmonicOscillator, s0, t, args=(k, m))
x1, x2, x3, x4, x5, x6, x7, x8, x9 = sol[:, 0], sol[:, 2], sol[:, 4], sol[:, 6], sol[:, 8], sol[:, 10], sol[:, 12], sol[:, 14], sol[:, 16]
q1 = np.sin(np.pi*1/10)*x1 + np.sin(np.pi*2/10)*x2 + np.sin(np.pi*3/10)*x3 + np.sin(np.pi*4/10)*x4 + np.sin(np.pi*5/10)*x5 + np.sin(np.pi*6/10)*x6 + np.sin(np.pi*7/10)*x7 + np.sin(np.pi*8/10)*x8 + np.sin(np.pi*9/10)*x9
q2 = np.sin(2*np.pi*1/10)*x1 + np.sin(2*np.pi*2/10)*x2 + np.sin(2*np.pi*3/10)*x3 + np.sin(2*np.pi*4/10)*x4 + np.sin(2*np.pi*5/10)*x5 + np.sin(2*np.pi*6/10)*x6 + np.sin(2*np.pi*7/10)*x7 + np.sin(2*np.pi*8/10)*x8 + np.sin(2*np.pi*9/10)*x9
q3 = np.sin(3*np.pi*1/10)*x1 + np.sin(3*np.pi*2/10)*x2 + np.sin(3*np.pi*3/10)*x3 + np.sin(3*np.pi*4/10)*x4 + np.sin(3*np.pi*5/10)*x5 + np.sin(3*np.pi*6/10)*x6 + np.sin(3*np.pi*7/10)*x7 + np.sin(3*np.pi*8/10)*x8 + np.sin(3*np.pi*9/10)*x9
q4 = np.sin(4*np.pi*1/10)*x1 + np.sin(4*np.pi*2/10)*x2 + np.sin(4*np.pi*3/10)*x3 + np.sin(4*np.pi*4/10)*x4 + np.sin(4*np.pi*5/10)*x5 + np.sin(4*np.pi*6/10)*x6 + np.sin(4*np.pi*7/10)*x7 + np.sin(4*np.pi*8/10)*x8 + np.sin(4*np.pi*9/10)*x9
q5 = np.sin(5*np.pi*1/10)*x1 + np.sin(5*np.pi*2/10)*x2 + np.sin(5*np.pi*3/10)*x3 + np.sin(5*np.pi*4/10)*x4 + np.sin(5*np.pi*5/10)*x5 + np.sin(5*np.pi*6/10)*x6 + np.sin(5*np.pi*7/10)*x7 + np.sin(5*np.pi*8/10)*x8 + np.sin(5*np.pi*9/10)*x9


fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', xlim=(0, L), ylim=(-L/5, 2*L/3))
ax.grid()

line1, = plt.plot([], [], 'ro-', animated=True)
line2, = plt.plot([], [], 'bo-', animated=True)
line3, = plt.plot([], [], 'go-', animated=True)
line4, = plt.plot([], [], 'yo-', animated=True)
line5, = plt.plot([], [], 'co-', animated=True)
line6, = plt.plot([], [], 'mo-', animated=True)

time_template = 'time = %.2fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
    return line1, line2, line3, line4, line5, line6, time_text
#    return line1, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line1.set_data([0, 
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
    line2.set_data([0, L/2 + q1[i]], [10, 10])
    line3.set_data([0, L/2 + q2[i]], [20, 20])
    line4.set_data([0, L/2 + q3[i]], [30, 30])
    line5.set_data([0, L/2 + q4[i]], [40, 40])
    line6.set_data([0, L/2 + q5[i]], [50, 50])
    time_text.set_text(time_template % (i*dt))
    return line1, line2, line3, line4, line5, line6, time_text
#    return line1, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 100/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

ani.save('./gif/damped_ladderdVibration.gif', writer='pillow', fps=fps)