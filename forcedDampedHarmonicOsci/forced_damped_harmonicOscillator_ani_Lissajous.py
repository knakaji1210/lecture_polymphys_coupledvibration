# ordinary differential equation of damped harmonic oscillator

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def forcedHarmonicOscillator(s, t, k, m, c, Fa, af):

    x, v = s                # s = (x, v)
    dsdt = [v, (-k*x)/m - c*v/m + (Fa/m)*np.sin(af*t)]    # 粘性抵抗、強制振動の項を追加
    return dsdt

def getNearestIndex2value(list,value):
    index = np.abs(np.array(list) -value).argsort()[0].tolist()
    return index

# variables
try:
    k = float(input('spring constant [N/m] (default=10.0): '))
except ValueError:
    k = 10.0               # [N/m] spring constant
try:
    m = float(input('mass [kg] (default=1.0): '))
except ValueError:
    m = 1.0                # [kg] mass
try:
    c = float(input('damping coefficient [kg/s] (default=0.5): '))    # [kg/s] damping coefficient
except ValueError:
    c = 0.5
l = 20                      # [m] equilibrium length
afreq0 = np.sqrt(k/m)       # natural　angular frequency
rho = c/(2*m)
tau = 1/rho
print('check "rho = {0:.1f} 1/s < af = {1:.1f} 1/s"'.format(rho,afreq0))
if np.abs(rho - afreq0) < 0.05: # critical damping (rho == afreq0はほぼ無理)
    pass
elif rho - afreq0 > 0: # over damping
    pass
elif rho - afreq0 < 0: # under damping
    afreq = np.sqrt(afreq0**2 - rho**2)
    period = 2*np.pi/afreq      # period [s] (T)

tmax = 10*period             # [s] duration time
dt = 0.05                   # [s] interval time

# forced oscillation
try:
    Fa = float(input('amplitude for forced oscillation [N] (default=10.0): '))
except ValueError:
    Fa = 10.0
try:
    af_min = float(input('minimum anglar frequency for forced oscillation [1/s] (default=1.0): '))
except ValueError:
    af_min = 1.0
try:
    af_max = float(input('maximum anglar frequency for forced oscillation [1/s] (default=5.0): '))
except ValueError:
    af_max = 5.0
try:
    num_freq = int(input('number of anglar frequency (default=5): '))
except ValueError:
    num_freq = 5

af_list = np.linspace(af_min, af_max, num_freq)

# initial condition
x0 = 0.0
v0 = 0.0
s0 = [x0, v0]               # initial condition

t = np.arange(0, tmax, dt)

xi = []             # 周波数掃引の全ての入力信号を格納
xo = []             # 周波数掃引の全ての出力信号を格納
xamp = []           # 各周波数での出力振幅の最大値を格納
xpha = []           # 各周波数での出力信号の位相を格納
one_list = np.ones(len(t)).tolist()     # アニメーション用（各周波数でのアニメーション期間中に同じ数値をずっと表示させるため）
af_ani = []         # アニメーション用の入力周波数の最大値を格納
xamp_ani = []       # アニメーション用の出力振幅の最大値を格納
xpha_ani = []       # アニメーション用の出力信号の位相を格納

for af in af_list:
    af_ani.extend([n*af for n in one_list])         # 入力周波数を格納（アニメーション用）
    xi_af = Fa*np.sin(af*t)                         # 入力信号
    xi_af_last = xi_af[int(0.8*len(xi_af)):]        # 後半部分を抽出（前半は過渡応答を含むから）
    xi.extend(xi_af)                                # 入力信号（アニメーション用）
    sol = odeint(forcedHarmonicOscillator, s0, t, args=(k,m,c,Fa,af))  # ODEの解を求めている
    xo_af = sol[:, 0]                               # [x]が出てくる
    xo_af_last = xo_af[int(0.8*len(xo_af)):]        # 後半部分を抽出（前半は過渡応答を含むから）
    xo_max = np.max(xo_af_last)                     # 後半部分の最大値（最大振幅と見做す）
    xamp.append(xo_max)                             # 最大振幅を格納
    xamp_ani.extend([n*xo_max for n in one_list])   # 最大振幅を格納（アニメーション用）
    xo.extend(xo_af)                                # 出力信号（アニメーション用）
    ind = getNearestIndex2value(xo_af_last,0)       # 出力信号が0になるindexを抽出
    if af < afreq:
        xo_pha = (180/np.pi)*np.arcsin(np.abs(xi_af_last[ind])/Fa)
    else:
        xo_pha = 180 - (180/np.pi)*np.arcsin(np.abs(xi_af_last[ind])/Fa)
    xpha.append(xo_pha)                             # 出力位相を格納
    xpha_ani.extend([n*xo_pha for n in one_list])   # 最大振幅を格納（アニメーション用）   

xamp_max = np.max(xamp)

fig = plt.figure()
ax = fig.add_subplot(111, xlim=(-1.5*Fa, 1.5*Fa), ylim=(-1.5*xamp_max, 1.5*xamp_max))
ax.grid()
ax.set_axisbelow(True)
ax.set_xlabel('$x$ position [m]')
peri1_template = r'$\omega_r$ = {0:.2f} s$^{{-1}}$,'.format(afreq)
ax.text(0.1, 0.8, peri1_template, transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

lissajous, = plt.plot([], [], 'r', animated=True)
# ここでは[],[]としているが、下で***.set_data([0, l + x[i]], [0, 0])で実際の値を入れている

amp_o_template = r'$A_{{out}}$ = %.1f m'
amp_o_text = ax.text(0.8, 0.9, '', transform=ax.transAxes)

pha_o_template = r'$\theta_{{out}}$ = %.1f$\degree$'
pha_o_text = ax.text(0.8, 0.8, '', transform=ax.transAxes)

peri2_template = r'$\omega_f$ = %.2f s$^{{-1}}$'
peri2_text = ax.text(0.1, 0.7, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入

time_template = '$t$ = %.2f s'
time_text = ax.text(0.1, 0.9, '', transform=ax.transAxes) # 図形の枠を基準にした位置にテキストが挿入
# また、ここでは''としているが、下で time_text.set_textで実際のテキストを入れている

def init():                 # FuncAnimationでinit_funcで呼び出す
    time_text.set_text('')
    return lissajous, time_text, peri2_text, amp_o_text, pha_o_text

def update(i):              # ここのiは下のframes=fに対応した引数になっている
    lissajous.set_data([xi[i-100:i]],[xo[i-100:i]])
    time_text.set_text(time_template % (i*dt))
    peri2_text.set_text(peri2_template % af_ani[i])
    amp_o_text.set_text(amp_o_template % xamp_ani[i])
    pha_o_text.set_text(pha_o_template % xpha_ani[i])
    return lissajous, time_text, peri2_text, amp_o_text, pha_o_text

'''
y_triの中の重要部分は
x_tri1 = np.linspace(a, b,100)
のとき
(xtri - a)/(b - a)
になる 
'''

f = np.arange(0, len(af_list)*len(t))
frame_int = 1000 * dt       # [ms] interval between frames
fps = 1000/frame_int        # frames per second

ani = FuncAnimation(fig, update, frames=f,
                    init_func=init, blit=True, interval=frame_int, repeat=True)

savefile = './gif/forced_harmonicOsci_ani_lissajous_(k={0:.1f},m={1:.1f},c={2:.1f}).gif'.format(k,m,c)
ani.save(savefile, writer='pillow', fps=fps)

plt.show()