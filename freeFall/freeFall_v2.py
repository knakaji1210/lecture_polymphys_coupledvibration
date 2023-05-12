import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def freeFall(s, t, r, m):

    g  = 9.80665            # [m/s^2] gravitational acceleration
    y, vel = s              # s = (y, vel)
    dsdt = [vel, (-r*vel-m*g)/m]
    return dsdt

r = 50                      # air resistance
m = 100                     # [kg] mass
y0 = 0                      # [m]   initial position
v0 = 20                     # [m/s] initial velocity
s0 = [y0, v0]               # initial condition
t = 7                       # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

sol = odeint(freeFall, s0, t, args=(r,m))

fig = plt.figure()
ax = fig.add_subplot(111, xlabel='time [s]', ylabel='position [m] & velocity [m/s]', xlim=(0, 7), ylim=(-100, 50))
ax.grid()

ax.plot(t, sol[:, 0], 'b', label='position')  # plot of y
ax.plot(t, sol[:, 1], 'g', label='velocity')  # plot of v
ax.legend(loc='best')

savefile = "./png/freeFall.png"
fig.savefig(savefile, dpi=300)

plt.show()