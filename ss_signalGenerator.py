#!/usr/bin/env python
from pickletools import read_stringnl_noescape_pair

__author__ = 'toopazo'

import argparse
# import os
# from obspy.core import read
#from obspy.core import UTCDateTime
import numpy as np
# import matplotlib.pyplot as plt
# import math
from obspy.core import read
from obspy.core import Trace, Stream, UTCDateTime
import ss_fft
import ss_plot
from random import uniform, randint


def generate_signal(sigtype):
    print(sigtype)

    signal_station = "R000"
    signal_starttime = "2015-01-01T00:00:00"
    # 2) Generate signal
    if sigtype == "sine":
        ######################################################
        # Tessampling_interval for a perfect sine wave
        ######################################################
        sampling_rate = 100.0      # sampling rate
        sampling_interval = 1.0/sampling_rate     # sampling interval
        signal_elapsed_secs = 5
        t_time = np.arange(0, signal_elapsed_secs, sampling_interval)     # time vector
        signal_frequency = 5     # frequency of the signal [Hz]
        signal_amplitude = (2**23)
        # Si no pongo ruido hay problemas de -inf en log10 (problemas de precision float)
        y_signal = signal_amplitude*np.sin(2*np.pi*signal_frequency*t_time) + randint(-20, 0)
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

    #3) Convert channel ch4
    data = np.int32(y_signal)
    channel = 'SHZ'
    # Fill header attributes
    stasampling_interval = {'network': 'UNK', 'station': signal_station, 'location': 'UNK',
                            'channel': channel, 'npsampling_interval': len(data), 'sampling_rate': sampling_rate,
                            'mseed': {'dataquality': 'D'},
                            'starttime': UTCDateTime(str(signal_starttime))}
    st = Stream([Trace(data=data, header=stasampling_interval)])
    # outfile_name = signal_starttime + "_" + signal_station + "_" + channel + ".MSEED"
    # outfile_name = outfile_name.replace(":", "-")
    outfile_name = sigtype + ".MSEED"
    st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
    #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
    st1 = read(outfile_name)
    arg = "[generate_signal] MSEED created: %s" % st1[0]
    print(arg)
    return outfile_name


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--sigtype', action='store', help='signal to generate', required=True)
    #parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('--outfile', action='store', help='name of output file')
    # parser.add_argument('--day_plot', action='store_true', help='day plot of the given file(s), normally same channel')
    args = parser.parse_args()

    resp_outfile_name = generate_signal(sigtype=args.sigtype)
    print(resp_outfile_name)
    ss_fft.do_rfft(resp_outfile_name, resp_outfile_name.replace(".MSEED", "_fft.png"), 24)
    ss_plot.plot_file(resp_outfile_name, resp_outfile_name.replace(".MSEED", ".png"), None, None)