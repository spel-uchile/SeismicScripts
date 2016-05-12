#!/usr/bin/Python
# -*- coding: utf-8 -*-

__author__ = 'toopazo'


import subprocess
import os
# from os import listdir
# from os.path import isfile, join
import ss_utilities
from obspy.core import read


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
                infile_i2 = os.path.abspath(xxx_files[i+1])
                print "[apply_to_folder] Processing infile_i %s .." % infile_i
                ApplyToFolder.apply_to_file(infile1=infile_i, infile2=infile_i2)
                print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i
            except IndexError:
                pass
        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile1, infile2):
        # 1) Make sure user inputs are correct
        infile1 = os.path.normcase(infile1)
        infile1 = os.path.normpath(infile1)
        infile1 = os.path.realpath(infile1)
        print(infile1)
        infile2 = os.path.normcase(infile2)
        infile2 = os.path.normpath(infile2)
        infile2 = os.path.realpath(infile2)
        print(infile2)

        # 2) Construct Stream object
        st = read(infile1)
        st += read(infile2)
        # st = st.sort()

        # 3) outfile name
        filename, file_extension = os.path.splitext(infile1)
        outfile = str(filename)
        outfile += '_rChi'
        outfile += ".MSEED"

        # 4) Save file
        st.write(filename=outfile, format='MSEED')


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
