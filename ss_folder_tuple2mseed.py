#!/usr/bin/Python
# -*- coding: utf-8 -*-

__author__ = 'toopazo'


import subprocess
import os
# from os import listdir
# from os.path import isfile, join
import ss_utilities
import ss_tuple2mseed


class ApplyToFolder():
    def __init__(self):
        pass

    @staticmethod
    def apply_to_folder(infolder):
        print "[apply_to_folder] infolder %s " % infolder
        print "**********************************************************"

        # 1) uncompress ".xxx.gz" files
        xxxgz_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=".tuple.gz")
        for path_i in xxxgz_files:
            gz_path_i = os.path.abspath(path_i)
            print "[apply_to_folder] Uncompressing gz_path_i %s .." % gz_path_i
            cmdline = "gzip -d %s" % gz_path_i
            subprocess.call(cmdline, shell=True)    # resp = str(subprocess.call(cmdline, shell=True))
            # arg = "[convert_slist2mseed] cmdline %s, resp %s" % (cmdline, resp)
            # print arg

        print "**********************************************************"

        # 2) get ".xxx" files and apply "apply_to_file"
        xxx_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=".tuple")
        for path_i in xxx_files:
            infile_i = os.path.abspath(path_i)
            print "[apply_to_folder] Processing infile_i %s .." % infile_i
            ApplyToFolder.apply_to_file(infile=infile_i)
            print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i

        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile):
        # ss_tuple2mseed.tuple2mseed(infile=infile,
        #                            user_ch1='NC',
        #                            user_ch2='SHE',
        #                            user_ch3='SHN',
        #                            user_ch4='SHZ')
        ss_tuple2mseed.tuple2mseed(infile=infile,
                                   user_ch1='SHZ',
                                   user_ch2='SHN',
                                   user_ch3='SHE',
                                   user_ch4='NC')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('directory', help='directory to use', action='store')
    args = parser.parse_args()

    ApplyToFolder.apply_to_folder(infolder=args.directory)
