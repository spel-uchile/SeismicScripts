#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from os import listdir
from os.path import isfile, join
import subprocess

parser = argparse.ArgumentParser(description='Obspy wrapper: Apply \"spectrogram\" operation for infolder')
# parser.add_argument('--infolder', action='store', help='files to process', required=True)
parser.add_argument('directory', help='directory to use', action='store')
parser.add_argument('--str1', action='store', help='str1 to filter', default='MSEED')    # , required=True)
parser.add_argument('--str2', action='store', help='str2 to filter', default='MSEED')    # , required=True)
# parser.add_argument('--showplot', action='store_true', help='show plot instead of saving it')
# parser.add_argument('--dayplot', action='store_true', help='dayplot of files')
parser.add_argument('--filter', action='store', help='filter signal before ploting')

args = parser.parse_args()

# 0) De-crease porcess priority: This should avoid new streamspqr and new report collisions
pid = os.getpid()
cmdline = "renice +30 -p %s" % pid
resp = str(subprocess.call(cmdline, shell=True))
arg = "renice operation: resp = %s" % resp
print arg

# 1) Make sure user inputs are correct
filter_plot = args.filter
# dayplot = args.dayplot
# showplot = args.showplot
str1 = args.str1
print(str1)
str2 = args.str2
print(str2)
# Convert to real (no symlink) and full path
infolder_path = args.directory
infolder_path = os.path.normcase(infolder_path)
infolder_path = os.path.normpath(infolder_path)
infolder_path = os.path.realpath(infolder_path)
print(infolder_path)
# Get all files in folder that contain "str1" and "str2"
onlyfiles = [f for f in listdir(infolder_path) if isfile(join(infolder_path, f))]
onlyfiles.sort()
sel_files = []
for file_i in onlyfiles:
    if (str1 in str(file_i)) and (str2 in str(file_i)):
        sel_files.append(file_i)
sel_files.sort()
infolder_files = sel_files
print(infolder_files)


# 3) Plot
outfile_extension = '.png'
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
#         #st.plot(type=plot_option_type, outfile=outfile_name, size=(800, 600))
#         tr0_header = str(st[0])[:15]
#         st.spectrogram(log=True, title=tr0_header + '    ' + str(st[0].stats.starttime),
#                        outfile=outfile_name)
#     else:
#         #st.plot(type=plot_option_type, method='full')
#         tr0_header = str(st[0])[:15]
#         st.spectrogram(log=True, title=tr0_header + '    ' + str(st[0].stats.starttime))

#else:
for file_i in infolder_files:
    # Construct Stream object, dayplotidually
    st = read(file_i)
    print(st[0].stats)
    print(st[0].data)

    if filter_plot is not None:
        st.filter("lowpass", freq=int(filter_plot), corners=10)   # , zerophase=True

    filename, file_extension = os.path.splitext(file_i)
    plot_option_type = 'normal'
    outfile_name = str(filename)
    outfile_name += '_spectrogram'
    outfile_name += outfile_extension
    #st.plot(type=plot_option_type, outfile=outfile_name, size=(800, 600))
    tr0_header = str(st[0])[:15]
    st.spectrogram(log=True, title=tr0_header + '    ' + str(st[0].stats.starttime),
                   outfile=outfile_name)
