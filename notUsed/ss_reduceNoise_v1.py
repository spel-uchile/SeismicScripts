#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
#import matplotlib
#matplotlib.rcParams.chunksize(100000)

parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--outfile', action='store', help='name of output file')
parser.add_argument('--dayplot', action='store_true', help='dayplot of the given file(s), normally same channel')
parser.add_argument('--filter', action='store', help='Filter signal before ploting')
args = parser.parse_args()

# 1) Make sure user inputs are correct
filter_plot = args.filter
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

# 2.1) Substract Nois v1
print st.traces
print st[0]
st.detrend()
noise_val = 100
# for i in range(0, len(st[0].data)):
#     if st[0].data[i] >= 0
#         st[0].data[i] += -noise_val
#     else:
#         st[0].data[i] += + noise_val
for i in range(0, len(st[0].data)):
    if 0 <= abs(st[0].data[i]) <= noise_val:
        st[0].data[i] = 0
print st[0].stats
print st[0].data
# exit(0)

# 3) Plot
if dayplot:
    plot_option_type = 'dayplot'
else:
    plot_option_type = 'normal'
if outfile_plot is not None:
    outfile_plot_name = str(outfile_plot)
    outfile_plot_name += '.png'
    if filter_plot is not None:
        st.filter("lowpass", freq=int(filter_plot), corners=10)   # , zerophase=True
    st.plot(type=plot_option_type, outfile=outfile_plot_name, size=(800, 600))
else:
    print(st[0].stats)
    print(st[0].data)
    if filter_plot is not None:
        st.filter("lowpass", freq=int(filter_plot), corners=10)   # , zerophase=True
    st.plot(type=plot_option_type, method='full')
