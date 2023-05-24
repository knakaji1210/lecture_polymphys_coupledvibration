import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def harmonicOscillator(s, t, k, m):

    x1, v1, x2, v2, x3, v3, x4, v4, x5, v5, x6, v6, x7, v7, x8, v8, x9, v9, x10, v10,  = s
    a1 = -2*k*x1/m + k*x2/m
    a2 = k*x1/m - 2*k*x2/m + k*x3/m
    a3 = k*x2/m - 2*k*x3/m + k*x4/m
    a4 = k*x3/m - 2*k*x4/m + k*x5/m
    a5 = k*x4/m - 2*k*x5/m + k*x6/m
    a6 = k*x5/m - 2*k*x6/m + k*x7/m
    a7 = k*x6/m - 2*k*x7/m + k*x8/m
    a8 = k*x7/m - 2*k*x8/m + k*x9/m
    a9 = k*x8/m - 2*k*x9/m + k*x10/m        
    a10 = k*x9/m - 2*k*x10/m
    dsdt = [v1, a1, v2, a2, v3, a3, v4, a4, v5, a5, v6, a6, v7, a7, v8, a8, v9, a9, v10, a10]
    return dsdt

k = 10                      # [N/m] spring constant
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
l11 = 10                    # [m] equilibrium length
L = l1 + l2 + l3 + l4 + l5 + l6 + l7 + l8 + l9 + l10 + l11
x1_0 = 10                   # [m]   initial position
x2_0 = 0                    # [m]   initial position
x3_0 = 0                    # [m]   initial position
x4_0 = 0                    # [m]   initial position
x5_0 = 0                    # [m]   initial position
x6_0 = 0                    # [m]   initial position
x7_0 = 0                    # [m]   initial position
x8_0 = 0                    # [m]   initial position
x9_0 = 0                    # [m]   initial position
x10_0 = 0                   # [m]   initial position
v1_0 = 0                    # [m/s] initial velocity
v2_0 = 0                    # [m/s] initial velocity
v3_0 = 0                    # [m/s] initial velocity
v4_0 = 0                    # [m/s] initial velocity
v5_0 = 0                    # [m/s] initial velocity
v6_0 = 0                    # [m/s] initial velocity
v7_0 = 0                    # [m/s] initial velocity
v8_0 = 0                    # [m/s] initial velocity
v9_0 = 0                    # [m/s] initial velocity
v10_0 = 0                    # [m/s] initial velocity
s0 = [x1_0, v1_0, x2_0, v2_0, x3_0, v3_0, x4_0, v4_0, x5_0, v5_0, x6_0, v6_0, x7_0, v7_0, x8_0, v8_0, x9_0, v9_0, x10_0, v10_0]
af = np.sqrt(2*k/m)         # angular frequency    
t = 8*np.pi/af             # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(harmonicOscillator, s0, t, args=(k, m))
x1, x2, x3, x4, x5, x6, x7, x8, x9, x10 = sol[:, 0], sol[:, 2], sol[:, 4], sol[:, 6], sol[:, 8], sol[:, 10], sol[:, 12], sol[:, 14], sol[:, 16], sol[:, 18]
q1 = (1/2)*(x1 + np.sqrt(2)*x2+x3)
q2 = (1/2)*(-np.sqrt(2)*x1 + np.sqrt(2)*x3)
q3 = (1/2)*(x1 - np.sqrt(2)*x2+x3)

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='x position [m]', xlim=(0, L), ylim=(-1, 4))
ax.grid()

line1, = plt.plot([], [], 'ro-', animated=True)
line2, = plt.plot([], [], 'bo-', animated=True)
line3, = plt.plot([], [], 'go-', animated=True)
line4, = plt.plot([], [], 'yo-', animated=True)

time_template = 'time = %.1fs'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes)

def init():
    time_text.set_text('')
#    return line1, line2, line3, line4, time_text
    return line1, time_text

def update(i):              # ここのiは下のframes=np.arange(0, len(t))に対応した引数になっている
    line1.set_data([0, 
            l1 + x1[i], 
            l1 + l2 + x2[i], 
            l1 + l2 + l3 + x3[i], 
            l1 + l2 + l3 +l4 + x4[i], 
            l1 + l2 + l3 +l4 + l5 + x5[i],
            l1 + l2 + l3 +l4 + l5 + l6 + x6[i], 
            l1 + l2 + l3 +l4 + l5 + l6 + l7+ x7[i], 
            l1 + l2 + l3 +l4 + l5 + l6 + l7 + l8 + x8[i],
            l1 + l2 + l3 +l4 + l5 + l6 + l7 + l8 + l9 + x9[i],
            l1 + l2 + l3 +l4 + l5 + l6 + l7 + l8 + l9 + l10 + x10[i],
             L], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
#    line2.set_data([0, L/2 + q1[i]], [1, 1])
#    line3.set_data([0, L/2 + q2[i]], [2, 2])
#    line4.set_data([0, L/2 + q3[i]], [3, 3])
    time_text.set_text(time_template % (i*dt))
#    return line1, line2, line3, line4, time_text
    return line1, time_text

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=np.arange(0, len(t)),
                    init_func=init, blit=True, interval=frame_int, repeat=True)
plt.show()

#ani.save('./gif/ladderdVibration.gif', writer='pillow', fps=fps)