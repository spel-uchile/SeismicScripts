#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
# import os
# from obspy.core import read
#from obspy.core import UTCDateTime
import numpy as np
# import math
from obspy.core import read
from obspy.core import Trace, Stream, UTCDateTime
import ss_fft
import ss_plot
from random import uniform, randint
import matplotlib.pyplot as plt
import math


class MyFunctions():
    def __init__(self):
        pass

    @staticmethod
    def galois_lfsr(seed, taps):
        seed = np.uint16(seed)
        taps = np.uint16(taps)
        start_state = np.uint16(seed)      # Any nonzero start state will work
        arg = "start_state \t\t%s" % bin(start_state)
        print(arg)
        lf_state_register = np.uint16(start_state)
        period = 0
        bit_sequence = []
        while True:
            lsb = lf_state_register & 1             # /* Get LSB (i.e., the output bit). */
            bit_sequence.append(lsb)
            lf_state_register >>= 1                 # /* Shift register */
            lf_state_register ^= (-lsb) & taps    # /* If the output bit is 1, apply toggle mask.
                                                    # * The value has 1 at bits corresponding
                                                    # * to taps, 0 elsewhere. */
            period += 1
            # arg = "lf_state_register \t%s" % bin(lf_state_register)
            # print(arg)
            if lf_state_register == start_state:
                break

        # arg = "bit_sequence[0] %s" % bit_sequence[0]
        # print(arg)
        # arg = "bit_sequence[1] %s" % bit_sequence[1]
        # print(arg)
        return bit_sequence


    @staticmethod
    def lfsr(seed, taps):
        state_register, xor = seed, 0
        bit_sequence = []
        while True:
            for t in taps:
                xor += int(state_register[t-1])
            if xor % 2 == 0.0:
                xor = 0
            else:
                xor = 1
            # print xor
            bit_sequence.append(xor)
            state_register, xor = str(xor) + state_register[:-1], 0
            # print state_register
            if state_register == seed:
                break
        return bit_sequence


def generate_signal(sigtype):
    print(sigtype)

    signal_station = "R000"
    signal_starttime = "2015-01-01T00:00:00"
    y_signal = None
    sampling_rate = None
    # 2) Generate signal
    if sigtype == "sine":
        ######################################################
        # Tessampling_interval for a perfect sine wave
        ######################################################
        # time array
        sampling_rate = 100.0      # sampling rate
        sampling_interval = 1.0/sampling_rate     # sampling interval
        signal_elapsed_secs = 5
        t_time = np.arange(0, signal_elapsed_secs, sampling_interval)     # time vector

        # signal array
        signal_frequency = 5     # frequency of the signal [Hz]
        signal_amplitude = (2**23)
        # Si no pongo ruido hay problemas de -inf en log10 (problemas de precision float)
        # y_signal = signal_amplitude*np.sin(2*np.pi*signal_frequency*t_time) + randint(-20, 0)
        y_signal = signal_amplitude*np.sin(2*np.pi*signal_frequency*t_time)
        arg = "len(y_signal) = %s | y_signal = %s" % (len(y_signal), y_signal)
        print(arg)

        # dataonly = y_signal
        # samp_rate = sampling_rate
        # plt.plot(t_time, dataonly, marker='o', markersize=4)
        # plt.show()
        # plt.clf()
        ######################################################
        # Tessampling_interval for a perfect sine wave
        ######################################################

    elif sigtype == "prbs":
        # https://en.wikipedia.org/wiki/Linear-feedback_shift_register#Some_polynomials_for_maximal_LFSRs
        # 16bits 6543210987654321
        seed = 0b1010110011100001   # anyone except 0 will work
        # 16bits 6543210987654321
        # taps = 0b1011010000000000   # Polynom = x^16 + x^14 + x^13 + x^11 + 1
        # 16bits 6543210987654321
        taps = 0b1101000000001000   # Polynom = x^16 + x^15 + x^13 + x^4 + 1 is a maximal-length polynomial LFSR

        # 8bits    87654321
        # seed = 0b11001001   # anyone except 0 will work
        # taps = 0b11001001
        y_signal = MyFunctions.galois_lfsr(seed, taps)       # 2^10 = 1023

        # Minimum duty-cycle should be 0.5*Fs = 50 Hz => 0.2ms (and sampling every 10ms)
        arg = "len(y_signal) %s" % len(y_signal)    # arg = "len(y_signal) %s, y_signal %s" % (len(y_signal), y_signal)
        print(arg)
        arg = "y_signal[0] %s" % y_signal[0]
        print(arg)
        arg = "y_signal[1] %s" % y_signal[1]
        print(arg)
        arg = "y_signal[2] %s" % y_signal[2]
        print(arg)
        arg = "y_signal[3] %s" % y_signal[3]
        print(arg)
        arg = "y_signal[4] %s" % y_signal[4]
        print(arg)

        y_signal2 = []
        for sample_i in range(0, len(y_signal)):
            if sample_i % 10 == 0:
                val = y_signal[sample_i]
                y_signal2.append(val)
            # y_signal2.append(val)
            # y_signal2.append(val)
            # y_signal2.append(val)
            # y_signal2.append(val)
            # y_signal2.append(val)
        y_signal = y_signal2

        arg = "len(y_signal) %s" % len(y_signal)    # arg = "len(y_signal) %s, y_signal %s" % (len(y_signal), y_signal)
        print(arg)
        arg = "y_signal[0] %s" % y_signal[0]
        print(arg)
        arg = "y_signal[1] %s" % y_signal[1]
        print(arg)
        arg = "y_signal[2] %s" % y_signal[2]
        print(arg)
        arg = "y_signal[3] %s" % y_signal[3]
        print(arg)
        arg = "y_signal[4] %s" % y_signal[4]
        print(arg)

        # time array
        sampling_rate = 100.0
        sampling_interval = 1.0/sampling_rate     # sampling interval
        signal_elapsed_secs = len(y_signal)*sampling_interval
        t_time = np.arange(0, signal_elapsed_secs, sampling_interval)     # time vector

    elif sigtype == "discreteDirac":
        # time array
        sampling_rate = 100.0      # sampling rate
        sampling_interval = 1.0/sampling_rate     # sampling interval
        signal_elapsed_secs = 5
        t_time = np.arange(0, signal_elapsed_secs, sampling_interval)     # time vector

        # signal array
        signal_amplitude = 10    # (2**24)
        signal_mean = 3
        y_signal = []
        for sample_i in range(0, len(t_time)):
            y_signal.append(signal_mean)
        middle_i = int(len(t_time)/2)
        arg = "middle_i = %s" % middle_i
        print(arg)
        y_signal[middle_i] = signal_amplitude*1.0
        arg = "len(y_signal) = %s | y_signal = %s" % (len(y_signal), y_signal)
        print(arg)

    elif sigtype == "dc":
        # time array
        sampling_rate = 100.0      # sampling rate
        sampling_interval = 1.0/sampling_rate     # sampling interval
        signal_elapsed_secs = 5
        t_time = np.arange(0, signal_elapsed_secs, sampling_interval)     # time vector

        # signal array
        signal_amplitude = (2**23)
        # signal_amplitude = 749623
        y_signal = []
        for sample_i in range(0, len(t_time)):
                y_signal.append(signal_amplitude)

    elif sigtype == "square":
        # time array
        sampling_rate = 100.0      # sampling rate
        sampling_interval = 1.0/sampling_rate     # sampling interval
        signal_elapsed_secs = 5
        t_time = np.arange(0, signal_elapsed_secs, sampling_interval)     # time vector

        # signal array
        signal_amplitude = (2**23)
        ud_bit = 1
        y_signal = []
        for sample_i in range(0, len(t_time)):
            if sample_i % 50 == 0:
                ud_bit *= -1
            y_signal.append(ud_bit*signal_amplitude)

    else:
        arg = "Unknown sigtype %s" % sigtype
        print(arg)
        exit(0)

    #3) Save to MSEED in channel SHZ
    data = np.int32(y_signal)
    channel = 'SHZ'
    # Fill header attributes
    stasampling_interval = {'network': 'UNK', 'station': signal_station, 'location': 'UNK',
                            'channel': channel, 'npsampling_interval': len(data), 'sampling_rate': sampling_rate,
                            'mseed': {'dataquality': 'D'},
                            'starttime': UTCDateTime(str(signal_starttime))}
    st = Stream([Trace(data=data, header=stasampling_interval)])
    # outfile = signal_starttime + "_" + signal_station + "_" + channel + ".MSEED"
    # outfile = outfile.replace(":", "-")
    outfile = sigtype + ".MSEED"
    st.write(outfile, format='MSEED', encoding=11, reclen=256, byteorder='>')
    #st.write(outfile, format='MSEED', encoding=0, reclen=256)
    st1 = read(outfile)
    arg = "[generate_signal] MSEED created: %s" % st1[0]
    print(arg)

    # # 4) Plot
    # # Plot Amp-Freq and Phase-Freq graphs
    # fig, ax = plt.subplots(3, 1)
    #
    # # Create arrays to plot Amp-Freq and Phase-Freq graphs
    # # Amplitude
    # AmpY_arr.append(AmpY)
    # log10_AmpY = 20*math.log10(AmpY)
    # log10_AmpY_arr.append(log10_AmpY)
    # # Frequency axis
    # f_arr.append(f)
    # log10_f = 10*math.log10(f)
    # log10_f_arr.append(log10_f)
    # # Pahse
    # PhaseY = cmath.phase(Y)
    # PhaseY = PhaseY*180.0/math.pi
    # PhaseY_arr.append(PhaseY)
    #
    # # Plot
    # ax[0].plot(log10_f_arr, log10_AmpY_arr, marker='o')
    # # ax[0].set_title('Frequency Response')
    # ax[0].set_xlabel('Frequency [log Hz]')
    # ax[0].set_ylabel('Amplitude Response [dB]')
    # # Extra info (max Freq, dBFS,  etc)
    # # y_axis = norm_abs_rfft_dbfs_SARA
    # # x_axis = freq_SARA
    # # y_axis_argmax = y_axis.argmax()
    # # y_axis_max = y_axis.max()
    # # print('y_axis_argmax = %s' % y_axis_argmax)
    # # print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    # # print('y_axis_max() = %s' % y_axis_max)
    # # arg_freqinfo = ' %s [Hz]\n %s [dbFS]' % (x_axis[y_axis_argmax], y_axis_max)
    # # ax[0].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.003, y_axis_max),
    # #                xytext=(x_axis[y_axis_argmax]+3, y_axis_max*1.3),
    # #                arrowprops=dict(facecolor='black', shrink=0.05),
    # #                )
    #
    # ax[1].plot(log10_f_arr, log10_AmpY_arr, marker='o')
    # # ax[0].set_title('Frequency Response')
    # ax[1].set_xlabel('Frequency [log Hz]')
    # ax[1].set_ylabel('Amplitude Response [dB]')
    # # Extra info (max Freq, dBFS,  etc)
    # # y_axis = norm_abs_rfft_dbfs_SARA
    # # x_axis = freq_SARA
    # # y_axis_argmax = y_axis.argmax()
    # # y_axis_max = y_axis.max()
    # # print('y_axis_argmax = %s' % y_axis_argmax)
    # # print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    # # print('y_axis_max() = %s' % y_axis_max)
    # # arg_freqinfo = ' %s [Hz]\n %s [dbFS]' % (x_axis[y_axis_argmax], y_axis_max)
    # # ax[0].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.003, y_axis_max),
    # #                xytext=(x_axis[y_axis_argmax]+3, y_axis_max*1.3),
    # #                arrowprops=dict(facecolor='black', shrink=0.05),
    # #                )
    #
    # ax[2].plot(log10_f_arr, PhaseY_arr, marker='o')
    # ax[2].set_xlabel('Frequency [log Hz]')
    # ax[2].set_ylabel('Phase Response[Hz]')
    # # Extra info (max Freq, dBFS,  etc)
    # # y_axis = norm_abs_rfft_dbfs_SPQR
    # # x_axis = freq_SPQR
    # # y_axis_argmax = y_axis.argmax()
    # # y_axis_max = y_axis.max()
    # # print('y_axis_argmax = %s' % y_axis_argmax)
    # # print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    # # print('y_axis_max() = %s' % y_axis_max)
    # # arg_freqinfo = ' %s [Hz]\n %s [dbFS]' % (x_axis[y_axis_argmax], y_axis_max)
    # # ax[1].annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.003, y_axis_max),
    # #                xytext=(x_axis[y_axis_argmax]+3, y_axis_max*1.3),
    # #                arrowprops=dict(facecolor='black', shrink=0.05),
    # #                )
    #
    # outfile = outfile.replace(".MSEED", "_TimeAmpPhase.png")
    # plt.savefig(outfile)
    # plt.show()
    # plt.clf()

    return outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--sigtype', action='store', help='signal to generate', required=True)
    #parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('--outfile', action='store', help='name of output file')
    # parser.add_argument('--day_plot', action='store_true', help='day plot of the given file(s), normally same channel')
    args = parser.parse_args()

    generate_signal(sigtype=args.sigtype)
    # resp_outfile = generate_signal(sigtype=args.sigtype)
    # print(resp_outfile)
    # ss_fft.do_rfft(resp_outfile, resp_outfile.replace(".MSEED", "_fft.png"), 24)
    # ss_plot.plot_file(resp_outfile, resp_outfile.replace(".MSEED", ".png"), None, None)