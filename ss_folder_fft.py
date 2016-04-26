#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from os import listdir
from os.path import isfile, join

import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Obspy wrapper: Apply \"FFT\" operation for infolder')
parser.add_argument('--infolder', action='store', help='files to process', required=True)
parser.add_argument('--str1', action='store', help='str2 to filter', required=True)
parser.add_argument('--str2', action='store', help='str2 to filter', required=True)
# parser.add_argument('--showplot', action='store_true', help='show plot instead of saving it')
# parser.add_argument('--dayplot', action='store_true', help='dayplot of files')
# parser.add_argument('--filter', action='store', help='filter signal before ploting')

args = parser.parse_args()

# 1) Make sure user inputs are correct
# filter_plot = args.filter
# dayplot = args.dayplot
# showplot = args.showplot
str1 = args.str1
print(str1)
str2 = args.str2
print(str2)
# Convert to real (no symlink) and full path
infolder_path = args.infolder
infolder_path = os.path.normcase(infolder_path)
infolder_path = os.path.normpath(infolder_path)
infolder_path = os.path.realpath(infolder_path)
print(infolder_path)
# Get all files in folder that contain "str1" and "str2"
onlyfiles = [f for f in listdir(infolder_path) if isfile(join(infolder_path, f))]
onlyfiles.sort()
sel_files = []
for file_i in onlyfiles:
    if (str1 in str(file_i)) and (str2 in str(file_i)):
        sel_files.append(file_i)
sel_files.sort()
infolder_files = sel_files
print(infolder_files)

outfile_extension = '.png'
for file_i in infolder_files:
    # 2) Construct Stream object
    st = read(file_i)
    print(st[0].stats)
    print(st[0].data)

    dataonly = st[0].data
    samp_rate = st[0].stats.sampling_rate

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
    plt.ylabel('Normalized |Y(f)| relative to Full Scale [dbFS], 24 bits')
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

    # 7) Save to file
    filename, file_extension = os.path.splitext(file_i)
    outfile_name = str(filename)
    outfile_name += '_fft'
    outfile_name += outfile_extension

    # tr0_header = str(st[0])[:15]
    # plt.title(tr0_header + '    ' + str(st[0].stats.starttime))
    # # plt.title(tr0_header + '    ' + str(st[0].stats.starttime) + '\n' + arg_freqinfo)
    # # plt.text(40, -40, arg_freqinfo)
    # arg_freqinfo = '%s [Hz]\n %s [dBFS]' % (freq[spec_abs_db_argmax], spec_abs_db_max)
    # plt.annotate(arg_freqinfo, xy=(freq[spec_abs_db_argmax]+1, spec_abs_db_max),
    #              xytext=(freq[spec_abs_db_argmax]+4, spec_abs_db_max-8),
    #              arrowprops=dict(facecolor='black', shrink=0.05),
    #              )
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Normalized response [dB]')
    plt.savefig(outfile_name)
    plt.clf()


    #################3
    #
    # if filter_plot is not None:
    #     st.filter("lowpass", freq=int(filter_plot), corners=10)   # , zerophase=True
    #
    # filename, file_extension = os.path.splitext(file_i)
    # plot_option_type = 'normal'
    # outfile_name = str(filename)
    # outfile_name += '_fft'
    # outfile_name += outfile_extension
    #
    # dataonly = st[0].data
    #
    # spec = np.fft.rfft(dataonly)
    # samp_rate = st[0].stats.sampling_rate
    # print samp_rate
    # freq = np.fft.rfftfreq(len(dataonly), d=1./samp_rate)
    # print spec.size
    # print freq.size
    # # spec = spec[:50]
    # # freq = freq[:50]
    # print spec.size
    # print freq.size
    #
    # spec_abs = abs(spec)
    # spec_abs_argmax = spec_abs.argmax()
    # spec_abs_max = spec_abs.max()
    # print('spec_abs.argmax() = %s' % spec_abs_argmax)
    # print('spec_abs[%s] = %s' % (spec_abs_argmax, spec_abs[spec_abs_argmax]))
    # print('spec_abs.max() = %s' % spec_abs_max)
    # arg = 'Max frequency component is freq[%s] = %s [Hz] = %s [mHz]' % (spec_abs_argmax, freq[spec_abs_argmax], 1000.0*freq[spec_abs_argmax])
    # print(arg)
    #
    # #plt.plot(freq, abs(spec))
    # markerline, stemlines, baseline = plt.stem(freq, spec_abs, '-.')
    # plt.setp(markerline, 'markerfacecolor', 'b')
    # plt.setp(baseline, 'color', 'r', 'linewidth', 2)
    #
    # tr0_header = str(st[0])[:15]
    # plt.title(tr0_header + '    ' + str(st[0].stats.starttime))
    # plt.xlabel('Frequency [Hz]\n '+arg)
    # plt.savefig(outfile_name)
    # plt.clf()
    # #plt.close()
