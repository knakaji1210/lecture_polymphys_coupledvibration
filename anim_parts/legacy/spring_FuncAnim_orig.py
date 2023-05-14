import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#t = 1           # [s] duration time
#dt = 0.1        # [s] interval time

k = 10                     # [N/m] spring constant
m = 50                      # [kg] mass
l = 20                      # [m] equilibrium length
afreq = np.sqrt(k/m)        # angular frequency    
t = 4*np.pi/afreq           # [s] duration time
dt = 0.1                    # [s] interval time
ds = 0.01
smax = 6

amp = 5

t = np.arange(0, t+dt, dt)
s = np.arange(0, smax*np.pi, ds)

print(len(t))
print(np.arange(0, len(t)))

fig = plt.figure()
ax = fig.add_subplot(111)

def update(i):     #ここでは利用してない
    ax.cla() # ax をクリア
    x = [t for t in t]
    y = [np.sin(t+i) for t in t]
    ax.plot(x, y)

def plot(i):
    ax.cla()
    ax.set_xlabel('x position [m]')
    ax.set_xlim(0,50)
    ax.set_ylim(-5,5)
    ax.grid()
    x = [(l + amp*np.cos(4*np.pi*i/100))*(np.sin((11/2)*s - np.pi/2)+(4/np.pi)*s + 1)/(smax*4+2) for s in s]
    y = [np.cos((11/2)*s - np.pi/2) for s in s]
    ax.plot(x,y)
    

frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, plot, frames=np.arange(0, len(t)), interval=frame_int, repeat=True)

plt.show()

ani.save('./gif/spring_orig.gif', writer='pillow', fps=fps)