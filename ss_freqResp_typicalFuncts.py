#!/usr/bin/env python
__author__ = 'toopazo'
 # # # -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import math
import cmath
import argparse


def freq_response(sigtype, order):
    order = int(order)
    arg = "[freq_response] sigtype %s .." % sigtype
    print(arg)
    arg = "[freq_response] order %s .." % order
    print(arg)
    y_amp = []
    y_phase = []
    frequencies = []
    freq_step = 0.0025
    sigtype_title = 'Empty Titel'
    freq_cutoff = 0
    max_freq = 44
    if sigtype == "my":
        for freq_i in range(1, int(max_freq/freq_step)):
            sigtype_title = 'My function \n H(jw) = '
            freq_val = freq_step*freq_i
            jw = 1j*(2.0*math.pi*freq_val)

            # F(s)
            freq_cutoff = 10
            w1 = 2.0*math.pi*freq_cutoff
            w2 = 2.0*math.pi*40
            # y_val = ((jw - 10)*(jw - 5))/((jw - 3))
            y_val = ((w1 / (jw + w1))**7)
            y_amp_val = abs(y_val)

            # Amplitude
            y_amp.append(y_amp_val)

            # Pahse
            y_phase_val = cmath.phase(y_val)
            y_phase_val = y_phase_val*180.0/math.pi
            y_phase.append(y_phase_val)

            # Frequency axis
            frequencies.append(freq_val)

            # Print info arround cut-off frequency
            # if freq_cutoff - 0.5 <= freq_val <= freq_cutoff + 0.5:
            if freq_cutoff == freq_val or freq_val == 1 or freq_val % 10 == 0:
                arg = "[freq_response] log10(%s) = %s [Hz] | 20*log10|%s| = %s [dB] | Phase(y) = %d"\
                      % (freq_val, np.log10(freq_val), y_amp_val, 20*np.log10(y_amp_val), y_phase_val)
                print(arg)

    elif sigtype == "regSARA":
        for freq_i in range(1, int(max_freq/freq_step)):
            sigtype_title = 'SARA Transfer Function\n H(jw)'
            freq_val = freq_step*freq_i
            jw = 1j*(2.0*math.pi*freq_val)

            # F(jw)
            freq_cutoff = 1.0

            y_val = (float('2.41083317014e-10')) * (jw - 563.732843616j) * (jw - 282.250600078j) * \
                    (jw - (-132.438457998+103.355980485j)) * (jw - (132.438457998+103.355980485j))
            # y_val *= (2**-23)
            # y_val *= (2**23)
            # y_val *= (2**24)/4.0
            # y_val /= (4.0)
            # y_val *= (2**-23)
            y_val *= (2**24)/4.0
            y_amp_val = abs(y_val)

            # Amplitude
            y_amp.append(y_amp_val)

            # Pahse
            y_val = 1.0     # there is no phase delay in ADC systems
            y_phase_val = cmath.phase(y_val)
            y_phase_val = y_phase_val*180.0/math.pi
            y_phase.append(y_phase_val)

            # Frequency axis
            frequencies.append(freq_val)

            # Print info arround cut-off frequency
            # if freq_cutoff - 0.5 <= freq_val <= freq_cutoff + 0.5:
            if freq_cutoff == freq_val or freq_val == 1 or freq_val % 10 == 0:
                arg = "[freq_response] log10(%s) = %s [Hz] | 20*log10|%s| = %s [dB] | Phase(y) = %d"\
                      % (freq_val, np.log10(freq_val), y_amp_val, 20*np.log10(y_amp_val), y_phase_val)
                print(arg)

    elif sigtype == "trillium40":
        for freq_i in range(1, int(max_freq/freq_step)):
            sigtype_title = 'Trillium 40 Transfer Function\n H(jw)'
            freq_val = freq_step*freq_i
            jw = 1j*(2.0*math.pi*freq_val)

            # F(jw)
            freq_cutoff = 1.0
            y_gain = float('1.104e5')   # * 1553.0  # gain 1553.0 [V / m/s]
            y_zeros = (jw - 0) * (jw - 0) * (jw + 68.8) * (jw + 323) * (jw + 2530)
            y_poles = (jw - (-0.1103 + 1j*0.1110)) * (jw - (-0.1103 - 1j*0.1110)) * \
                      (jw + 86.3) * \
                      (jw - (-241 + 1j*178)) * (jw - (-241 - 1j*178)) * \
                      (jw - (-535 + 1j*719)) * (jw - (-535 - 1j*719))
            y_val = y_gain*y_zeros/y_poles
            y_amp_val = abs(y_val)

            # Amplitude
            y_amp.append(y_amp_val)

            # Pahse
            y_phase_val = cmath.phase(y_val)
            y_phase_val = y_phase_val*180.0/math.pi
            y_phase.append(y_phase_val)

            # Frequency axis
            frequencies.append(freq_val)

            # Print info arround cut-off frequency
            # if freq_cutoff - 0.5 <= freq_val <= freq_cutoff + 0.5:
            if freq_cutoff == freq_val or freq_val == 1 or freq_val % 10 == 0:
                arg = "[freq_response] log10(%s) = %s [Hz] | 20*log10|%s| = %s [dB] | Phase(y) = %d"\
                      % (freq_val, np.log10(freq_val), y_amp_val, 20*np.log10(y_amp_val), y_phase_val)
                print(arg)

    elif sigtype == "lp":
        for freq_i in range(1, int(max_freq/freq_step)):
            sigtype_title = 'LP filter \n H(jw) = (wo / (jw + wo))**%s' % order
            freq_val = freq_step*freq_i
            jw = 1j*(2.0*math.pi*freq_val)

            # F(s) of a LP filter, see also https://stanford.edu/~boyd/ee102/conv_demo.pdf
            freq_cutoff = 10
            wo = 2.0*math.pi*freq_cutoff
            y_val = (wo / (jw + wo))**order
            y_amp_val = abs(y_val)

            # Amplitude
            y_amp.append(y_amp_val)

            # Pahse
            y_phase_val = cmath.phase(y_val)
            y_phase_val = y_phase_val*180.0/math.pi
            y_phase.append(y_phase_val)

            # Frequency axis
            frequencies.append(freq_val)

            # Print info arround cut-off frequency
            # if freq_cutoff - 0.5 <= freq_val <= freq_cutoff + 0.5:
            if freq_cutoff == freq_val or freq_val == 1 or freq_val % 10 == 0:
                arg = "[freq_response] log10(%s) = %s [Hz] | 20*log10|%s| = %s [dB] | Phase(y) = %d"\
                      % (freq_val, np.log10(freq_val), y_amp_val, 20*np.log10(y_amp_val), y_phase_val)
                print(arg)

    elif sigtype == "hp":
        for freq_i in range(1, int(max_freq/freq_step)):
            sigtype_title = 'HP filter \n H(jw) = (jw / (jw + wo))**%s' % order
            freq_val = freq_step*freq_i
            jw = 1j*(2.0*math.pi*freq_val)

            # F(s) of a HP filter, see also https://stanford.edu/~boyd/ee102/conv_demo.pdf
            freq_cutoff = 10
            wo = 2.0*math.pi*freq_cutoff
            y_val = (jw / (jw + wo))**order
            y_amp_val = abs(y_val)

            # Amplitude
            y_amp.append(y_amp_val)

            # Pahse
            y_phase_val = cmath.phase(y_val)
            y_phase_val = y_phase_val*180.0/math.pi
            y_phase.append(y_phase_val)

            # Frequency axis
            frequencies.append(freq_val)

            # Print info arround cut-off frequency
            # if freq_cutoff - 0.5 <= freq_val <= freq_cutoff + 0.5:
            if freq_cutoff == freq_val or freq_val == 1 or freq_val % 10 == 0:
                arg = "[freq_response] log10(%s) = %s [Hz] | 20*log10|%s| = %s [dB] | Phase(y) = %d"\
                      % (freq_val, np.log10(freq_val), y_amp_val, 20*np.log10(y_amp_val), y_phase_val)
                print(arg)

    elif sigtype == "dampedPendulum":
        for freq_i in range(1, int(max_freq/freq_step)):
            sigtype_title = 'Damped Pendulum \n H(jw) = (jw)*(jw) / (jw-p1)*(jw-p2)'
            freq_val = freq_step*freq_i
            jw = 1j*(2.0*math.pi*freq_val)

            # F(s) of a damped pendulum
            d = 0.707
            freq_cutoff = 10.0
            p1 = d + 1j*math.sqrt((1-d**2)**2)
            p1 = p1*2.0*math.pi*freq_cutoff
            p2 = d - 1j*math.sqrt((1-d**2)**2)
            p2 = p2*2.0*math.pi*freq_cutoff
            y_val = ((jw - 0)*(jw - 0))/((jw-p1)*(jw-p2))
            y_amp_val = abs(y_val)

            # Amplitude
            y_amp.append(y_amp_val)

            # Pahse
            y_phase_val = cmath.phase(y_val)
            y_phase_val = y_phase_val*180.0/math.pi
            y_phase.append(y_phase_val)

            # Frequency axis
            frequencies.append(freq_val)

            # Print info arround cut-off frequency
            # if freq_cutoff - 0.5 <= freq_val <= freq_cutoff + 0.5:
            if freq_cutoff == freq_val or freq_val == 1 or freq_val % 10 == 0:
                arg = "[freq_response] log10(%s) = %s [Hz] | 20*log10|%s| = %s [dB] | Phase(y) = %d"\
                      % (freq_val, np.log10(freq_val), y_amp_val, 20*np.log10(y_amp_val), y_phase_val)
                print(arg)

    # Create arrays to plot Amp-Freq and Phase-Freq graphs
    # Amplitude
    y_amp_db = 20*np.log10(y_amp)
    # Frequency axis
    frequencies_log10 = 10*np.log10(frequencies)

    # Plot Amp-Freq and Phase-Freq graphs
    plt.subplot(2, 1, 1)
    # print len(ax)

    # ax[0].plot(frequencies_log10, y_amp_db, marker='o')
    plt.plot(frequencies, y_amp_db, marker='o')
    plt.title(sigtype_title)
    # plt.xlabel('Frequency [Hz]')
    # ax[0].xscale('log', basey=10, nonposy='clip')
    # plt.xscale('log', basey=10, nonposy='clip')
    plt.ylabel('Amplitude Response [dB]')
    # Extra info (max Freq, dBFS,  etc)
    y_axis = np.float64(y_amp_db)
    x_axis = np.float64(frequencies)
    y_axis_argmax = int(freq_cutoff/freq_step) - 1
    y_axis_max = y_axis[y_axis_argmax]
    y_axis_min = y_axis.min()
    x_axis_min = x_axis.min()
    print('y_axis_argmax = %s' % y_axis_argmax)
    print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    print('y_axis_max() = %s' % y_axis_max)
    arg_freqinfo = '%s [Hz]\n %s [dB]' % (x_axis[y_axis_argmax], y_axis_max)
    plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax], y_axis_max),
                 xytext=(x_axis_min+1, y_axis_min+1),
                 arrowprops=dict(facecolor='black', shrink=0.05),
                 )

    plt.subplot(2, 1, 2)
    # ax[1].plot(frequencies_log10, y_phase, marker='o')
    plt.plot(frequencies, y_phase, marker='o')
    plt.xlabel('Frequency [Hz]')
    # ax[1].xscale('log', basey=10, nonposy='clip')
    # plt.xscale('log', basey=10, nonposy='clip')
    plt.ylabel('Phase Response [degrees]')
    # Extra info (max Freq, dBFS,  etc)
    y_axis = y_phase
    x_axis = frequencies
    y_axis_argmax = int(freq_cutoff/freq_step) - 1
    y_axis_max = y_axis[y_axis_argmax]
    print('y_axis_argmax = %s' % y_axis_argmax)
    print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    print('y_axis_max() = %s' % y_axis_max)
    arg_freqinfo = ' %s [Hz]\n %s [degrees]' % (x_axis[y_axis_argmax], y_axis_max)
    plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+1, y_axis_max+1),
                 xytext=(x_axis[y_axis_argmax]+5, y_axis_max+10),
                 arrowprops=dict(facecolor='black', shrink=0.05),
                 )

    if order != 1:
        outfile = "trasferFunction_" + sigtype + "_order" + str(order) + "_AmpPhase.png"
    else:
        outfile = "trasferFunction_" + sigtype + "_AmpPhase.png"
    arg = "[freq_response] saving file %s .." % outfile
    print(arg)
    plt.savefig(outfile)
    # plt.show()
    plt.clf()

    return frequencies, y_amp_db, y_phase

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--sigtype', action='store', help='signal to generate', required=True)
    parser.add_argument('--order', action='store', help='function order')
    #parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('--outfile', action='store', help='name of output file')
    # parser.add_argument('--day_plot', action='store_true', help='day plot of the given file(s), normally same channel')
    args = parser.parse_args()

    if args.order is None:
        args.order = 1
    freq_response(sigtype=args.sigtype, order=args.order)