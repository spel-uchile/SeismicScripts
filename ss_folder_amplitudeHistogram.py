#!/usr/bin/Python
# -*- coding: utf-8 -*-

__author__ = 'toopazo'


import subprocess
import os
# from os import listdir
# from os.path import isfile, join
import numpy as np
from obspy.core import read, Stream, Trace
import matplotlib.pyplot as plt
# import plotly.plotly as py

import ss_utilities


class ApplyToFolder():
    def __init__(st):
        pass

    @staticmethod
    def apply_to_folder(infolder, bins_width, bins_max):
        #>print "[apply_to_folder] infolder %s " % infolder
        #>print "**********************************************************"

        # 1) uncompress ".xxx.gz" files
        xxxgz_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=".MSEED.gz")
        for path_i in xxxgz_files:
            gz_path_i = os.path.abspath(path_i)
            #>print "[apply_to_folder] Uncompressing gz_path_i %s .." % gz_path_i
            cmdline = "gzip -d %s" % gz_path_i
            subprocess.call(cmdline, shell=True)    # resp = str(subprocess.call(cmdline, shell=True))
            # arg = "[convert_slist2mseed] cmdline %s, resp %s" % (cmdline, resp)
            # #>print arg

        #>print "**********************************************************"

        # 2) get ".xxx" files and apply "apply_to_file"
        xxx_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=".MSEED")
        for i in range(0, len(xxx_files)):
            try:
                infile_i = os.path.abspath(xxx_files[i])
                #>print "[apply_to_folder] Processing infile_i %s .." % infile_i
                ApplyToFolder.apply_to_file(infile=infile_i, bins_width=bins_width, bins_max=bins_max)
                #>print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i
            except IndexError:
                pass
        #>print "**********************************************************"

        # 3) Done
        #>print "Done"

    @staticmethod
    def apply_to_file(infile, bins_width, bins_max):
        # 1) Make sure user inputs are correct
        infile = os.path.normcase(infile)
        infile = os.path.normpath(infile)
        infile = os.path.realpath(infile)
        #>print(infile)
        bins_width = int(bins_width)
        bins_max = int(bins_max)

        # 2) Construct Stream object
        st = read(infile)
        # st = st.sort()

        # tr_i) Histogram
        bins_array = []
        bins_i = -bins_max + bins_width/2
        while bins_i <= bins_max:
            bins_array.append(bins_i)
            bins_i += bins_width
        # #>print(bins_array)
        for tr_i in range(0, len(st)):
            arg = 'st[tr_i].stats.channel %s' % st[tr_i].stats.channel
            # #>print(arg)
            print(st[tr_i].stats)
            samples = st[tr_i].data
            samples = np.int32(samples)

            # #>print(st[tr_i])
            # #>print(st[tr_i].data)

            # Histogram example from https://plot.ly/matplotlib/histograms/
            freq_of_bins_normed, bins, patches = plt.hist(samples, bins=bins_array, normed=True)
            # arg = "len(freq_of_bins_normed) %s freq_of_bins_normed %s"\
            #       % (len(freq_of_bins_normed), freq_of_bins_normed)
            # #>print(arg)
            # arg = "len(bins) %s bins %s" % (len(bins), bins)
            # #>print(arg)
            # arg = "len(patches) %s patches %s" % (len(patches), patches)
            # #>print(arg)

            # get statistical info about sampels
            samples_mean = np.mean(samples)
            arg = 'samples_mean %s' % samples_mean
            print(arg)
            samples_std = np.std(samples)
            arg = 'samples_std %s' % samples_std
            print(arg)
            samples_max = np.max(samples)
            arg = 'samples_max %s' % samples_max
            print(arg)
            samples_min = np.min(samples)
            arg = 'samples_min %s' % samples_min
            print(arg)
            pseudo_signal = 0
            pseudo_noise = 0
            for sample_i in range(0, len(samples)):
                if abs(samples[sample_i]) <= 3.0*samples_std:
                    pseudo_noise += abs(samples[sample_i])
                else:
                    pseudo_signal += abs(samples[sample_i])
            try:
                pseudo_snr = 1.0*pseudo_signal/pseudo_noise
            except ZeroDivisionError:
                pseudo_snr = "ZeroDivisionError"
            arg = 'pseudo_snr %s' % pseudo_snr
            print(arg)
            arg = "Histogram \n Mean %s stdDev %s \n pseudo_snr %s" % (samples_mean, samples_std, pseudo_snr)
            plt.title(arg)

            # get statistical info about bins
            plt.yscale('log', basey=10, nonposy='clip')
            plt.xlabel("Amplitude intervals [normal scale]")
            plt.ylabel("Monrmed Frequency [log10 scale]")
            x_axis = np.float64(bins_array)
            y_axis = np.float64(freq_of_bins_normed)
            y_axis_argmax = y_axis.argmax()
            y_axis_max = y_axis.max()
            #>print('y_axis_argmax = %s' % y_axis_argmax)
            #>print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
            #>print('y_axis_max = %s' % y_axis_max)
            # arg_freqinfo = 'Maximum: \n%s [Interval]\n %s [Frequency]' % (x_axis[y_axis_argmax], y_axis_max)
            # plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+0.0, y_axis_max+0.0),
            #              xytext=(0.5, 0.5),
            #              arrowprops=dict(facecolor='black', shrink=0.05),
            #              textcoords='figure fraction'
            #              )
            # plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+10, y_axis_max+0.0),
            #              xytext=(x_axis[y_axis_argmax]+100, y_axis_max),
            #              #arrowprops=dict(facecolor='black', shrink=0.05),
            #              )

            # outfile_channel_ext = "_" + st[tr_i].stats.channel + "_ampHisto.png"
            outfile = infile.replace(".MSEED", "_ampHisto.png")
            plt.savefig(outfile)
            # plt.show()
            plt.clf()
            print("")
        # exit(0)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('directory', help='directory to use', action='store')
    parser.add_argument('--bins_width', action='store', help='bins_width', required=True)
    parser.add_argument('--bins_max', action='store', help='bins_max', required=True)
    args = parser.parse_args()

    uinfolder = args.directory
    uinfolder = os.path.normcase(uinfolder)
    uinfolder = os.path.normpath(uinfolder)
    uinfolder = os.path.realpath(uinfolder)

    ApplyToFolder.apply_to_folder(uinfolder, bins_width=args.bins_width, bins_max=args.bins_max)
