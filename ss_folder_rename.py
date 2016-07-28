#!/usr/bin/env python
__author__ = 'toopazo'

import os
# import argparse
# from obspy.core import read
# from os import listdir
# from os.path import isfile, join
import ss_plot
import ss_utilities
import subprocess


class ApplyToFolder():
    def __init__(self):
        pass

    @staticmethod
    def apply_to_folder(infolder, ext, strin, strout):
        print "[apply_to_folder] infolder %s " % infolder
        print "**********************************************************"

        # 1) get ".xxx" files and apply "apply_to_file"
        xxx_files = ss_utilities.ParserUtilities.get_xxx_files(folderpath=infolder, extension=ext)
        for path_i in xxx_files:
            infile_i = os.path.abspath(path_i)
            print "[apply_to_folder] Processing infile_i %s .." % infile_i
            ApplyToFolder.apply_to_file(infile_i, strin, strout)
            print "Next file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n"    # separate infile_i

        print "**********************************************************"

        # 3) Done
        print "Done"

    @staticmethod
    def apply_to_file(infile, strin, strout):
        # basename = os.path.basename(infile)
        # print "[apply_to_file] basename %s .." % basename
        directory, filename = os.path.split(infile)
        print "[apply_to_file] directory %s filename %s .." % (directory, filename)

        filename = filename.replace(strin, strout)
        outfile = directory + "/" + filename
        outfile = os.path.normcase(outfile)
        outfile = os.path.normpath(outfile)
        outfile = os.path.realpath(outfile)
        print "[apply_to_file] filename %s .." % filename

        os.rename(infile, outfile)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Apply function xxx to corresponding files .yyy in the given folder')
    parser.add_argument('directory', help='directory to use', action='store')
    parser.add_argument('--strin', action='store', help='string in', required=True)     # , default='MSEED')
    parser.add_argument('--strout', action='store', help='string out', required=True)     # , default='MSEED')
    parser.add_argument('--ext', action='store', help='extension to filter', default='MSEED')
    args = parser.parse_args()

    uinfolder = args.directory
    uinfolder = os.path.normcase(uinfolder)
    uinfolder = os.path.normpath(uinfolder)
    uinfolder = os.path.realpath(uinfolder)

    ApplyToFolder.apply_to_folder(infolder=uinfolder, ext=args.ext, strin=args.strin, strout=args.strout)


# parser = argparse.ArgumentParser(description='Obspy wrapper: Apply \"plot\" operation for infolder')
# parser.add_argument('--infolder', action='store', help='files to process', required=True)
# parser.add_argument('--str1', action='store', help='str2 to filter', required=True)
# parser.add_argument('--str2', action='store', help='str2 to filter', required=True)
# parser.add_argument('--showplot', action='store_true', help='show plot instead of saving it')
# parser.add_argument('--dayplot', action='store_true', help='dayplot of files')
# parser.add_argument('--filter', action='store', help='filter signal before ploting')
#
# args = parser.parse_args()
#
# # 1) Make sure user inputs are correct
# filter_plot = args.filter
# dayplot = args.dayplot
# showplot = args.showplot
# str1 = args.str1
# print(str1)
# str2 = args.str2
# print(str2)
# # Convert to real (no symlink) and full path
# infolder_path = args.infolder
# infolder_path = os.path.normcase(infolder_path)
# infolder_path = os.path.normpath(infolder_path)
# infolder_path = os.path.realpath(infolder_path)
# print(infolder_path)
# # Get all files in folder that contain "str1" and "str2"
# onlyfiles = [f for f in listdir(infolder_path) if isfile(join(infolder_path, f))]
# onlyfiles.sort()
# sel_files = []
# for file_i in onlyfiles:
#     if (str1 in str(file_i)) and (str2 in str(file_i)):
#         sel_files.append(file_i)
# sel_files.sort()
# infolder_files = sel_files
# print(infolder_files)
#
#
# # 3) Plot
# outfile_extension = '.png'
# if dayplot:
#     # Construct Stream object, appending every trace in the folder
#     st = read(infolder_files[0])
#     for i in range(1, len(infolder_files)):
#         st += read(infolder_files[i])
#     st = st.sort()
#     print(st[0].stats)
#     print(st[0].data)
#
#     if filter_plot is not None:
#         st.filter("lowpass", freq=int(filter_plot), corners=10)   # , zerophase=True
#
#     plot_option_type = 'dayplot'
#     if showplot is not None:
#         outfile_name = 'dayplot'
#         outfile_name += outfile_extension
#         st.plot(type=plot_option_type, outfile=outfile_name, size=(800, 600))
#     else:
#         st.plot(type=plot_option_type, method='full')
#
# else:
#     for file_i in infolder_files:
#         # Construct Stream object, dayplotidually
#         st = read(file_i)
#         print(st[0].stats)
#         print(st[0].data)
#
#         if filter_plot is not None:
#             st.filter("lowpass", freq=int(filter_plot), corners=10)   # , zerophase=True
#
#         filename, file_extension = os.path.splitext(file_i)
#         plot_option_type = 'normal'
#         outfile_name = str(filename)
#         outfile_name += outfile_extension
#         st.plot(type=plot_option_type, outfile=outfile_name, size=(800, 600))
