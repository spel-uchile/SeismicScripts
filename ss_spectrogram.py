#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
#from obspy.core import UTCDateTime

parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--outfile', action='store', help='name of output file')
# parser.add_argument('--day_plot', action='store_true', help='day plot of the given file(s), normally same channel')
args = parser.parse_args()

# 1) Make sure user inputs are correct
outfile_plot = args.outfile
outfile_plot_name = str(outfile_plot)
outfile_plot_name += '.png'
# day_plot = args.day_plot
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
# if day_plot:
#     plot_option_type = 'dayplot'
# else:
#     plot_option_type = 'normal'
if outfile_plot is not None:
    outfile_plot_name = str(outfile_plot)
    outfile_plot_name += '.png'
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(800, 600))
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(1024, 768))
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(1920, 1080))
    #st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(2048, 1080))
    #st.filter("lowpass", freq=7, corners=10)   # , zerophase=True
    tr0_header = str(st[0])[:15]
    st.spectrogram(log=True, title=tr0_header + '    ' + str(st[0].stats.starttime),
                   outfile=outfile_plot_name)
else:
    #st.plot(type=plot_option_type, method='full')
    tr0_header = str(st[0])[:15]
    # st.filter("lowpass", freq=7, corners=10)   # , zerophase=True
    st.spectrogram(log=True, title=tr0_header + '    ' + str(st[0].stats.starttime),
                   outfile=outfile_plot_name)
