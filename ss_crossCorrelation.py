#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mplt
# mplt.rcParams.chunksize(100000)
# mplt.rcParams['chunksize'] = (100000)
# mplt.rcParams['agg.path.chunksize'] = 10000
# mplt.use('Agg')


def cross_correlation(infile1, infile2, outfile):
    # 1) Make sure user inputs are correct
    # Convert to real (no symlink) and full path
    infile1 = os.path.normcase(infile1)
    infile1 = os.path.normpath(infile1)
    infile1 = os.path.realpath(infile1)
    print(infile1)
    infile2 = os.path.normcase(infile2)
    infile2 = os.path.normpath(infile2)
    infile2 = os.path.realpath(infile2)
    print(infile1)
    if outfile is not None:
        outfile = os.path.normcase(outfile)
        outfile = os.path.normpath(outfile)
        outfile = os.path.realpath(outfile)
    else:
        # ss_utilities.ParserUtilities.
        infile1 = infile1.split("/")
        infile1 = infile1[len(infile1)-1]
        print(infile1)
        infile2 = infile2.split("/")
        infile2 = infile2[len(infile2)-1]
        print(infile2)
        outfile = infile1.replace(".MSEED", "_and_") + infile2.replace(".MSEED", "_crossCorr.png")
    print(outfile)

    # 2) Construct Stream object
    st1 = read(infile1)
    samples1 = np.float64(st1[0].data)
    # Normalization according to https://en.wikipedia.org/wiki/Cross-correlation#Normalized_cross-correlation
    samples1 -= np.mean(samples1)
    samples1 /= np.std(samples1)
    samples1 /= len(samples1)
    arg = "len(samples1) %s samples1 %s" % (len(samples1), samples1)
    print(arg)

    st2 = read(infile2)
    samples2 = np.float64(st2[0].data)
    samples2 -= np.mean(samples2)
    samples2 /= np.std(samples2)
    arg = "len(samples2) %s samples2 %s" % (len(samples2), samples2)
    print(arg)

    cross_corr = np.correlate(samples1, samples2, mode="full")
    # cross_corr_max = cross_corr.max()
    # cross_corr /= cross_corr_max    # Normalize relative to the maximum point
    # cross_corr /= cross_corr[len(cross_corr)/2]    # Normalize relative to the center point
    arg = "len(cross_corr) %s cross_corr %s" % (len(cross_corr), cross_corr)
    print(arg)

    # time array
    # sampling_rate = 100.0      # sampling rate
    # sampling_interval = 1.0/sampling_rate     # sampling interval
    sampling_interval = 1
    starting_index = -len(cross_corr)/2 + 1
    ending_index = len(cross_corr)/2 + 1
    t_time = np.arange(starting_index, ending_index, sampling_interval)     # time vector
    arg = "len(t_time) %s t_time %s" % (len(t_time), t_time)
    print(arg)

    # Plot
    # cross_corr = np.log10(cross_corr)
    plt.plot(t_time, cross_corr, marker='o', color='black', linestyle="-", markersize=4)
    # plt.yscale('log', basey=10, nonposy='clip')

    # # Extra info (max Freq, dBFS,  etc)
    plt.title('Cross-Correlation')
    plt.xlabel('Samples Lag [Sampling Period]')
    plt.ylabel('Normalized Cross-Correlation')
    # y_axis = cross_corr
    # x_axis = t_time
    # y_axis_argmax = y_axis.argmax()
    # y_axis_max = y_axis.max()
    # print('y_axis_argmax = %s' % y_axis_argmax)
    # print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    # print('y_axis_max() = %s' % y_axis_max)
    # arg_freqinfo = ' %s [Sample Lag]\n %s [Cross-Correlation]' % (x_axis[y_axis_argmax], y_axis_max)
    # plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.1, y_axis_max),
    #              xytext=(x_axis[y_axis_argmax]+100, y_axis_max-0.10),
    #              arrowprops=dict(facecolor='black', shrink=0.05),
    #              )

    # 7) Save to file
    plt.savefig(outfile)
    # plt.show()
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile1', action='store', help='files to process', required=True)
    parser.add_argument('--infile2', action='store', help='files to process', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file')
    args = parser.parse_args()

    cross_correlation(infile1=args.infile1,
                      infile2=args.infile2,
                      outfile=args.outfile)
