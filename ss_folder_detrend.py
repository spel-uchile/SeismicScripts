#!/usr/bin/Python
# -*- coding: utf-8 -*-

__author__ = 'toopazo'


import subprocess
import os
# from os import listdir
from obspy import read
from obspy.core import read, Stream, Trace
import numpy as np

import ss_utilities
import ss_spqrSampleCorrection


class ApplyToFolder():
    def __init__(self):
        pass

    @staticmethod
    def apply_to_folder(infolder, udc, utrend):
        print "[apply_to_folder] infolder %s " % infolder
        print "**********************************************************"

        # 1) uncompress ".xxx.gz" files
        xxxgz_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=".MSEED.gz")
        for path_i in xxxgz_files:
            gz_path_i = os.path.abspath(path_i)
            print "[apply_to_folder] Uncompressing gz_path_i %s .." % gz_path_i
            cmdline = "gzip -d %s" % gz_path_i
            subprocess.call(cmdline, shell=True)    # resp = str(subprocess.call(cmdline, shell=True))
            # arg = "[convert_slist2mseed] cmdline %s, resp %s" % (cmdline, resp)
            # print arg

        print "**********************************************************"

        # 2) get ".xxx" files and apply "apply_to_file"
        xxx_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=".MSEED")
        for path_i in xxx_files:
            infile_i = os.path.abspath(path_i)
            print "[apply_to_folder] Processing infile_i %s .." % infile_i
            ApplyToFolder.apply_to_file(infile=infile_i, dc=udc, trend=utrend)
            print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i

        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile, dc, trend):
        # 1) Make sure user inputs are correct
        infile = os.path.normcase(infile)
        infile = os.path.normpath(infile)
        infile = os.path.realpath(infile)
        arg = "[apply_to_file] infile %s" % infile
        print(arg)
        arg = "[apply_to_file] dc %s" % dc
        print(arg)
        arg = "[apply_to_file] trend %s" % trend
        print(arg)

        # 2) Construct Stream object
        st = read(infile)

        if trend:
            # 4) Filter and detrend
            for tr_i in range(0, len(st)):
                st[tr_i].data = st[tr_i].detrend().data
                # st[tr_i].filter("lowpass", freq=40, corners=10)
                pass

        if dc:
            # 3) Substract average (self-made DC detrend, actually works better)
            for tr_i in range(0, len(st)):
                samples = st[tr_i].data
                # samples = np.int32(samples)
                samples_mean = np.mean(samples)
                arg = "[apply_to_file] pre samples_mean %s" % samples_mean
                print(arg)
                # samples_mean = samples.mean()
                samples -= samples_mean
                # samples *= 0
                samples_mean = np.mean(samples)
                arg = "[apply_to_file] post samples_mean %s" % samples_mean
                print(arg)
                # for j in range(0, len(samples)):
                #     samples[j] -= samples_mean
                st[tr_i].data = samples
                # print(st[tr_i])
                # print(st[tr_i].data)

        # last step) Save file
        st.write(filename=infile, format='MSEED')
        # st = Stream([Trace(data=st.data, header=st.stats)])
        # st.write(infile, format='MSEED', encoding=11, reclen=256, byteorder='>')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('--dc', action='store_true', help='substract DC', default=True)
    parser.add_argument('--trend', action='store_true', help='substract trend')
    # parser.add_argument('--str1', action='store', help='str1 to filter', default='MSEED')    # , required=True)
    # parser.add_argument('--str2', action='store', help='str2 to filter', default='MSEED')    # , required=True)
    parser.add_argument('directory', help='directory to use', action='store')
    args = parser.parse_args()

    uinfolder = args.directory
    uinfolder = os.path.normcase(uinfolder)
    uinfolder = os.path.normpath(uinfolder)
    uinfolder = os.path.realpath(uinfolder)
    udc = args.dc
    utrend = args.trend

    ApplyToFolder.apply_to_folder(uinfolder, udc, utrend)
