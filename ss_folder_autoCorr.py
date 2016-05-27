#!/usr/bin/env python
__author__ = 'toopazo'

import os
# import argparse
# from obspy.core import read
# from os import listdir
# from os.path import isfile, join
import ss_crossCorrelation
import ss_utilities
import subprocess


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
        for path_i in xxx_files:
            infile_i = os.path.abspath(path_i)
            print "[apply_to_folder] Processing infile_i %s .." % infile_i
            ApplyToFolder.apply_to_file(infile_i, str1, str2)
            print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i

        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile, str1, str2):
        if (str1 in str(infile)) and (str2 in str(infile)):
            outfile = infile.replace(".MSEED", "_autoCorr.png")
            ss_crossCorrelation.cross_correlation(infile1=infile, infile2=infile, outfile=outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('directory', help='directory to use', action='store')
    parser.add_argument('--str1', action='store', help='str2 to filter', required=True)
    parser.add_argument('--str2', action='store', help='str2 to filter', required=True)
    args = parser.parse_args()

    uinfolder = args.directory
    uinfolder = os.path.normcase(uinfolder)
    uinfolder = os.path.normpath(uinfolder)
    uinfolder = os.path.realpath(uinfolder)

    ApplyToFolder.apply_to_folder(infolder=uinfolder,
                                  str1=args.str1,
                                  str2=args.str2)
