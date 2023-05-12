import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

t = 7                       # [s] duration time
dt = 0.1                    # [s] interval time

t = np.arange(0, t+dt, dt)

fig = plt.figure()
ax = fig.add_subplot(111)

def update(i):
    ax.cla() # ax をクリア
    x = [t for t in t]
    y = [np.sin(t+i) for t in t]
    ax.plot(x, y)

anim = FuncAnimation(fig, update, frames=np.arange(0, len(t)), interval=1000)

plt.show()

plt.close()