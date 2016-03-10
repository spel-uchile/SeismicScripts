#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
#from obspy.core import UTCDateTime
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--outfile', action='store', help='name of output file')
# parser.add_argument('--day_plot', action='store_true', help='day plot of the given file(s), normally same channel')
args = parser.parse_args()

# 1) Make sure user inputs are correct
outfile_plot = args.outfile
# day_plot = args.day_plot
# Convert to real (no symlink) and full path
infile_paths = args.infile
for i in range(0, len(infile_paths)):
    infile_paths[i] = os.path.normcase(infile_paths[i])
    infile_paths[i] = os.path.normpath(infile_paths[i])
    infile_paths[i] = os.path.realpath(infile_paths[i])
    print(infile_paths[i])

# 2) Construct Stream object
st = read(infile_paths[0])
for i in range(1, len(infile_paths)):
    st += read(infile_paths[i])
st = st.sort()
print(st[0].stats)
print(st[0].data)

# 3) Calculate FFT
dataonly = st[0].data
samp_rate = st[0].stats.sampling_rate
# ######################################################
# # Tests for a perfect sine wave
# ######################################################
# import math
# Fs = 100.0      # sampling rate
# Ts = 1.0/Fs     # sampling interval
# t_time = np.arange(0, 1, Ts)     # time vector
# ff = 5     # frequency of the signal [Hz]
# y_signal = (2**23)*np.sin(2*np.pi*ff*t_time)
# dataonly = y_signal
# samp_rate = Fs
# # plt.plot(t_time, dataonly, marker='o', markersize=4)
# # plt.show()
# # plt.clf()
# ######################################################
# # Tests for a perfect sine wave
# ######################################################
spec = np.fft.rfft(dataonly)
freq = np.fft.rfftfreq(len(dataonly), d=1./samp_rate)
print('samp_rate %s' % samp_rate)
print('spec %s' % spec)
print('freq %s' % freq)
print('spec.size %s' % spec.size)
print('freq.size %s' % freq.size)

# 4) Generate Normalized Fourier transform it to dBFS
norm_abs_rfft = 2.0*abs(spec)/len(dataonly)
norm_abs_rfft_dbfs = norm_abs_rfft
len_norm_abs_rfft = len(norm_abs_rfft)
nbits = 24
for i in range(0, len_norm_abs_rfft):
    ratio = (1.0*norm_abs_rfft[i])/(2.0**(nbits-1))
    # ratio = (1.0*norm_abs_rfft[i])
    norm_abs_rfft_dbfs[i] = 20.0*np.log10(ratio)
print('norm_abs_rfft_dbfs (24-bit ADC is assumed !!) %s' % norm_abs_rfft_dbfs)

# 5) Plot
plt.plot(freq, norm_abs_rfft_dbfs, marker='o', markersize=4)
plt.xlabel('Freq (Hz)')
plt.ylabel('Normalized |Y(f)| relative to Full Scale [dbFS]')
# Extra info (max Freq, dBFS,  etc)
y_axis = norm_abs_rfft_dbfs
x_axis = freq
y_axis_argmax = y_axis.argmax()
x_axis_max = y_axis.max()
print('y_axis_argmax = %s' % y_axis_argmax)
print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
print('x_axis_max() = %s' % x_axis_max)
arg_freqinfo = ' %s [Hz]\n %s [dBFS]' % (x_axis[y_axis_argmax], x_axis_max)
plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+1, x_axis_max),
             xytext=(x_axis[y_axis_argmax]+4, x_axis_max-10),
             arrowprops=dict(facecolor='black', shrink=0.05),
             )

# 6) Save or show
if outfile_plot is not None:
    outfile_plot_name = str(outfile_plot)
    outfile_plot_name += '.png'
    plt.savefig(outfile_plot_name)
else:
    plt.show()




# ######################################################
# # Tests for a perfect sine wave
# ######################################################
# Fs = 100.0      # sampling rate
# Ts = 1.0/Fs     # sampling interval
# t_time = np.arange(0, 1, Ts)     # time vector
# ff = 5     # frequency of the signal [Hz]
# y_signal = (2**23)*np.sin(2*np.pi*ff*t_time)
#
# fig, ax = plt.subplots(3, 1)
# ax[0].plot(t_time, y_signal, marker='o')
# ax[0].set_xlabel('Time')
# ax[0].set_ylabel('y(t)')
# # Extra info (max Freq, dBFS,  etc)
# y_axis = y_signal
# x_axis = t_time
# y_axis_argmax = y_axis.argmax()
# x_axis_max = y_axis.max()
# print('y_axis_argmax = %s' % y_axis_argmax)
# print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
# print('x_axis_max() = %s' % x_axis_max)
# arg_freqinfo = ' %s [s]\n %s [Counts]' % (x_axis[y_axis_argmax], x_axis_max)
# ax[0].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.003, x_axis_max),
#                xytext=(x_axis[y_axis_argmax]+0.03, x_axis_max*0.7),
#                arrowprops=dict(facecolor='black', shrink=0.05),
#                )
#
# spec = np.fft.rfft(y_signal)
# freq = np.fft.rfftfreq(len(y_signal), d=Ts)
# norm_fft = spec/len(y_signal)
# norm_abs_rfft = 2.0*abs(norm_fft)
# ax[1].plot(freq, norm_abs_rfft, marker='o', markersize=4)
# ax[1].set_xlabel('Freq (Hz)')
# ax[1].set_ylabel('Normalized |Y(f)|')
# # Extra info (max Freq, dBFS,  etc)
# y_axis = norm_abs_rfft
# x_axis = freq
# y_axis_argmax = y_axis.argmax()
# x_axis_max = y_axis.max()
# print('y_axis_argmax = %s' % y_axis_argmax)
# print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
# print('x_axis_max() = %s' % x_axis_max)
# arg_freqinfo = ' %s [Hz]\n %s [Counts]' % (x_axis[y_axis_argmax], x_axis_max)
# ax[1].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+1, x_axis_max),
#                xytext=(x_axis[y_axis_argmax]+4, x_axis_max*0.7),
#                arrowprops=dict(facecolor='black', shrink=0.05),
#                )
#
# # Generate dBFS
# norm_abs_rfft_dbfs = norm_abs_rfft
# len_norm_abs_rfft = len(norm_abs_rfft)
# nbits = 24
# for i in range(0, len_norm_abs_rfft):
#     # norm = spec_abs[i]*(1.0/len_norm_abs_rfft)
#     if norm_abs_rfft[i] >= (1000):
#         print('norm_abs_rfft[%s] = %s !!!!!!!!!' % (i, norm_abs_rfft[i]))
#     ratio = (1.0*norm_abs_rfft[i])/(2.0**(nbits-1))
#     # ratio = (1.0*spec_abs[i])
#     norm_abs_rfft_dbfs[i] = 20.0*np.log10(ratio)
# ax[2].plot(freq, norm_abs_rfft_dbfs, marker='o', markersize=4)
# ax[2].set_xlabel('Freq (Hz)')
# ax[2].set_ylabel('Normalized |Y(f)| [dbFS]')
# # Extra info (max Freq, dBFS,  etc)
# y_axis = norm_abs_rfft_dbfs
# x_axis = freq
# y_axis_argmax = y_axis.argmax()
# x_axis_max = y_axis.max()
# print('y_axis_argmax = %s' % y_axis_argmax)
# print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
# print('x_axis_max() = %s' % x_axis_max)
# arg_freqinfo = ' %s [Hz]\n %s [dBFS]' % (x_axis[y_axis_argmax], x_axis_max)
# ax[2].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+1, x_axis_max),
#                xytext=(x_axis[y_axis_argmax]+4, x_axis_max-70),
#                arrowprops=dict(facecolor='black', shrink=0.05),
#                )
#
# # Plot and erase
# plt.show()
# plt.clf()
#
# exit(0)
# ######################################################
# # Tests for a perfect sine wave
# ######################################################