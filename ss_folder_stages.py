#!/usr/bin/Python
# -*- coding: utf-8 -*-

__author__ = 'toopazo'


import subprocess
import os
# from os import listdir
# from os.path import isfile, join
import numpy as np
from obspy.core import read, Stream, Trace
import random

import ss_utilities


class ApplyToFolder():
    def __init__(self):
        pass

    @staticmethod
    def apply_to_folder(infolder, str1, str2):
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
        for i in range(0, len(xxx_files)):
            try:
                infile_i = os.path.abspath(xxx_files[i])
                print "[apply_to_folder] Processing infile_i %s .." % infile_i
                ApplyToFolder.apply_to_file(infile=infile_i, str1=str1, str2=str2)
                print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i
            except IndexError:
                pass
        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile, str1, str2):
        if (str1 in str(infile)) and (str2 in str(infile)):
            # 1) Make sure user inputs are correct
            infile = os.path.normcase(infile)
            infile = os.path.normpath(infile)
            infile = os.path.realpath(infile)
            print(infile)

            # >> Construct Stream object
            st = read(infile)

            for trace_i in range(0, len(st)):
                # >> Substract DC component
                st[trace_i].data = np.int32(st[trace_i].data - np.mean(st[trace_i].data))

                # >> lowpass filter
                # cutoff_freq = 45
                # st.filter("lowpass", freq=int(cutoff_freq), corners=10)   # , zerophase=True

                # >>Add AWGN noise
                # noise_std_dev = 10
                # arg = "[apply_to_file] noise_std_dev %s" % noise_std_dev
                # print(arg)
                # samples = st[trace_i].data
                # samples = np.int32(samples)
                # for sample_i in range(0, len(samples)):
                #     samples[sample_i] += random.gauss(0, noise_std_dev)
                # st[trace_i].data = samples

                # >>

            # >> Save changes
            st.write(infile, format='MSEED', encoding=11, reclen=256, byteorder='>')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('--str1', action='store', help='str1 to filter', default='MSEED')    # , required=True)
    parser.add_argument('--str2', action='store', help='str2 to filter', default='MSEED')    # , required=True)
    parser.add_argument('directory', help='directory to use', action='store')
    args = parser.parse_args()

    uinfolder = args.directory
    uinfolder = os.path.normcase(uinfolder)
    uinfolder = os.path.normpath(uinfolder)
    uinfolder = os.path.realpath(uinfolder)

    ApplyToFolder.apply_to_folder(infolder=uinfolder, str1=args.str1, str2=args.str2)

