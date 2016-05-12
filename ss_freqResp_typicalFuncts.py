#!/usr/bin/env python
__author__ = 'toopazo'

import numpy as np
import matplotlib.pyplot as plt
import math
import cmath

log10_AmpY_arr = []
AmpY_arr = []
PhaseY_arr = []
f_arr = []
log10_f_arr = []
for i in range(1, 100*3):
    f = 1.0*i/3
    jw = 1j*(2.0*math.pi*f)

    # F(s) of a damped pendulum
    d = 0.707
    fo = 10.0
    p1 = d + 1j*math.sqrt((1-d**2)**2)
    p1 = p1*2.0*math.pi*fo
    p2 = d - 1j*math.sqrt((1-d**2)**2)
    p2 = p2*2.0*math.pi*fo
    Y = ((jw**2)/((jw-p1)*(jw-p2)))
    AmpY = abs(Y)

    # F(s) of a LP filter
    # fo = 20
    # wo = 2.0*math.pi*fo
    # Y = 1/(1 + jw/(wo))
    # AmpY = abs(Y)

    # F(s) of a HP filter
    # fo = 20
    # wo = 2.0*math.pi*fo
    # Y = (jw/wo)/(1 + jw/wo)
    # AmpY = abs(Y)

    # Create arrays to plot Amp-Freq and Phase-Freq graphs
    # Amplitude
    AmpY_arr.append(AmpY)
    log10_AmpY = 20*math.log10(AmpY)
    log10_AmpY_arr.append(log10_AmpY)
    # Frequency axis
    f_arr.append(f)
    log10_f = 10*math.log10(f)
    log10_f_arr.append(log10_f)
    # Pahse
    PhaseY = cmath.phase(Y)
    PhaseY = PhaseY*180.0/math.pi
    PhaseY_arr.append(PhaseY)

    # Print info arround cut-off frequency
    if fo - 1.0 <= f <= fo + 1.0:
        arg = "f %s [Hz], log10(f) %s => log10|AmpY| %s [dB] and Phase(Y) %d" % (f, log10_f, log10_AmpY, PhaseY)
        print(arg)

# Plot Amp-Freq and Phase-Freq graphs
fig, ax = plt.subplots(2, 1)

ax[0].plot(log10_f_arr, log10_AmpY_arr, marker='o')
ax[0].set_title('Frequency Response')
ax[0].set_xlabel('Frequency [log Hz]')
ax[0].set_ylabel('Amplitude Response [dB]')
# Extra info (max Freq, dBFS,  etc)
# y_axis = norm_abs_rfft_dbfs_SARA
# x_axis = freq_SARA
# y_axis_argmax = y_axis.argmax()
# x_axis_max = y_axis.max()
# print('y_axis_argmax = %s' % y_axis_argmax)
# print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
# print('x_axis_max() = %s' % x_axis_max)
# arg_freqinfo = ' %s [Hz]\n %s [dbFS]' % (x_axis[y_axis_argmax], x_axis_max)
# ax[0].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.003, x_axis_max),
#                xytext=(x_axis[y_axis_argmax]+3, x_axis_max*1.3),
#                arrowprops=dict(facecolor='black', shrink=0.05),
#                )

ax[1].plot(log10_f_arr, PhaseY_arr, marker='o')
ax[1].set_xlabel('Frequency [log Hz]')
ax[1].set_ylabel('Phase Response[Hz]')
# Extra info (max Freq, dBFS,  etc)
# y_axis = norm_abs_rfft_dbfs_SPQR
# x_axis = freq_SPQR
# y_axis_argmax = y_axis.argmax()
# x_axis_max = y_axis.max()
# print('y_axis_argmax = %s' % y_axis_argmax)
# print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
# print('x_axis_max() = %s' % x_axis_max)
# arg_freqinfo = ' %s [Hz]\n %s [dbFS]' % (x_axis[y_axis_argmax], x_axis_max)
# ax[1].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.003, x_axis_max),
#                xytext=(x_axis[y_axis_argmax]+3, x_axis_max*1.3),
#                arrowprops=dict(facecolor='black', shrink=0.05),
#                )

plt.show()