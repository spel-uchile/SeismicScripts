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
        ss_tuple2mseed.tuple2mseed(infile=infile,
                                   user_ch1='NC',
                                   user_ch2='SHE',
                                   user_ch3='SHN',
                                   user_ch4='SHZ')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('directory', help='directory to use', action='store')
    args = parser.parse_args()

    ApplyToFolder.apply_to_folder(infolder=args.directory)


# class PostProcConverter():
#     def __init__(self):
#         pass
#
#     @staticmethod
#     def convert_folder_tuple2mseed(folderpath):
#         print "[convert_folder_tuple2mseed] folderpath %s " % folderpath
#         print "**********************************************************"
#
#         # 1) uncompress ".gz" files
#         gz_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=folderpath, extension=".tuple.gz")
#         for path_i in gz_files:
#             gz_path_i = os.path.abspath(path_i)
#             print "[convert_folder_tuple2mseed] Uncompressing gz_path_i %s .." % gz_path_i
#             cmdline = "gzip -d %s" % gz_path_i
#             subprocess.call(cmdline, shell=True)    # resp = str(subprocess.call(cmdline, shell=True))
#             # arg = "[convert_slist2mseed] cmdline %s, resp %s" % (cmdline, resp)
#             # print arg
#
#         print "**********************************************************"
#
#         # 2) get .tuple files and convert them to ".slist", and then to ".MSEED"
#         tuple_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=folderpath, extension=".tuple")
#         for path_i in tuple_files:
#             tuple_path_i = os.path.abspath(path_i)
#             print "[convert_folder_tuple2mseed] Processing tuple_path_i %s .." % tuple_path_i
#             PostProcConverter.convert_tuple2mseed(tuple_path=tuple_path_i)
#             print "Next \".tuple\" file .."    # separate tuple_path_i
#
#         print "**********************************************************"
#
#         # 3) Done
#         print "Done"
#
#     @staticmethod
#     def convert_tuple2mseed(tuple_path):
#         ss_tuple2mseed.tuple2mseed(infile_paths=[tuple_path],
#                                    user_ch1='NC',
#                                    user_ch2='SHE',
#                                    user_ch3='SHN',
#                                    user_ch4='SHZ')
#
#
# if __name__ == '__main__':
#     import argparse
#     parser = argparse.ArgumentParser(description='Convert \".tuple\" files into \".slist\" and later to \".mseed\"')
#
#     parser.add_argument('directory', help='directory to use', action='store')
#     args = parser.parse_args()
#
#     infolder = os.path.abspath(args.directory)
#     arg = "infolder %s" % infolder
#     print(arg)
#
#     PostProcConverter.convert_folder_tuple2mseed(infolder)
