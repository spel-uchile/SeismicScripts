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
    def apply_to_folder(infolder):
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
                ApplyToFolder.apply_to_file(infile=infile_i)
                print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i
            except IndexError:
                pass
        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile):
        # 1) Make sure user inputs are correct
        infile = os.path.normcase(infile)
        infile = os.path.normpath(infile)
        infile = os.path.realpath(infile)
        print(infile)

        # Add noise
        for noise_i in range(1, 10):

            # Construct Stream object
            st = read(infile)

            # Add gaussian noise
            noise_std_dev = noise_i*4
            arg = "[apply_to_file] noise_std_dev %s" % noise_std_dev
            print(arg)
            for i in range(0, len(st)):
                samples = st[i].data
                samples = np.int32(samples)
                for sample_i in range(0, len(samples)):
                    samples[sample_i] += random.gauss(0, noise_std_dev)
                st[i].data = samples
                # print(st[i])
                # print(st[i].data)

            # Save file
            infile_ext = "_%s.MSEED" % noise_std_dev
            outfile = infile
            outfile = outfile.replace(".MSEED", infile_ext)
            st.write(filename=outfile, format='MSEED')
            # st = Stream([Trace(data=st.data, header=st.stats)])
            # st.write(infile, format='MSEED', encoding=11, reclen=256, byteorder='>')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('directory', help='directory to use', action='store')
    args = parser.parse_args()

    uinfolder = args.directory
    uinfolder = os.path.normcase(uinfolder)
    uinfolder = os.path.normpath(uinfolder)
    uinfolder = os.path.realpath(uinfolder)

    ApplyToFolder.apply_to_folder(uinfolder)

