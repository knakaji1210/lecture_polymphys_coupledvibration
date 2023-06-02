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
elif rho - afreq0 > 0:      # over damping
    pass
elif rho - afreq0 < 0:      # under damping
    afreq = np.sqrt(afreq0**2 - rho**2)
    period = 2*np.pi/afreq  # period [s] (T)

tmax = 10*period            # [s] duration time
dt = 0.01                   # [s] interval time

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
    num_freq = int(input('number of anglar frequency (default=100): '))
except ValueError:
    num_freq = 100

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
ax1 = fig.add_subplot(111, xlim=(af_min, af_max), ylim=(0, 1.2*xamp_max))
ax1.set_axisbelow(True)
ax1.set_xlabel('$\omega_f$ /s$^{{-1}}$')
ax1.set_ylabel(r'$A_{{out}}$ /m')

var1_template = r'$k$ = {0:.1f} N/m'.format(k)
var2_template = r'$m$ = {0:.1f} kg'.format(m)
var3_template = r'$c$ = {0:.1f} kg/s'.format(c)
ax1.text(0.1, 0.95, var1_template, transform=ax1.transAxes) # 図形の枠を基準にした位置にテキストが挿入
ax1.text(0.1, 0.85, var2_template, transform=ax1.transAxes) # 図形の枠を基準にした位置にテキストが挿入
ax1.text(0.1, 0.75, var3_template, transform=ax1.transAxes) # 図形の枠を基準にした位置にテキストが挿入
peri1_template = r'$\tau$ = {0:.2f} s'.format(tau)
peri2_template = r'$\omega_r$ = {0:.2f} s$^{{-1}}$,'.format(afreq)
ax1.text(0.1, 0.65, peri1_template, transform=ax1.transAxes) # 図形の枠を基準にした位置にテキストが挿入
ax1.text(0.1, 0.55, peri2_template, transform=ax1.transAxes) # 図形の枠を基準にした位置にテキストが挿入

ax1.plot(af_list,xamp, 'r')
ax1.vlines(afreq, 0, 1.2*xamp_max, 'g', ls='dashed', lw=1)

ax2 = ax1.twinx()
ax2.grid(ls='dotted')
ax2.set_ylim(0,180)
ax2.set_ylabel(r'$\theta_{{out}}$ /$\degree$')
ax2.plot(af_list,xpha, 'b')
ax2.hlines(90, 0, 1.2*xamp_max, 'g', ls='dashed', lw=1)

savefile = './png/forced_harmonicOsci_freqSweep_(k={0:.1f},m={1:.1f},c={2:.1f}).png'.format(k,m,c)
fig.savefig(savefile, dpi=300)

plt.show()