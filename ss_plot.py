#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read

parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--outfile', action='store', help='name of output file')
parser.add_argument('--dayplot', action='store_true', help='dayplot of the given file(s), normally same channel')
args = parser.parse_args()

# 1) Make sure user inputs are correct
outfile_plot = args.outfile
dayplot = args.dayplot
# Convert to real (no symlink) and full path
infile_paths = args.infile
for i in range(0, len(infile_paths)):
    infile_paths[i] = os.path.normcase(infile_paths[i])
    infile_paths[i] = os.path.normpath(infile_paths[i])
    infile_paths[i] = os.path.realpath(infile_paths[i])
    print(infile_paths[i])

# 2) Construct Stream object
st = read(infile_paths[0])
for i in range(1, len(infile_paths)):
    st += read(infile_paths[i])
st = st.sort()

# 3) Plot
if dayplot:
    plot_option_type = 'dayplot'
else:
    plot_option_type = 'normal'
if outfile_plot is not None:
    outfile_plot_name = str(outfile_plot)
    outfile_plot_name += '.png'
    st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(800, 600))
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(1024, 768))
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(1920, 1080))
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(2048, 1080))
else:
    print(st[0].stats)
    print(st[0].data)
    #st.filter("lowpass", freq=49, corners=10)   # , zerophase=True
    st.plot(type=plot_option_type, method='full')
