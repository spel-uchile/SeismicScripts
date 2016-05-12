#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
#import matplotlib
#matplotlib.rcParams.chunksize(100000)


def plot_file(infile, outfile, outfilter, outdayplot):
    # 1) Make sure user inputs are correct
    # Convert to real (no symlink) and full path
    infile = os.path.normcase(infile)
    infile = os.path.normpath(infile)
    infile = os.path.realpath(infile)
    print(infile)
    if outfile is not None:
        outfile = os.path.normcase(outfile)
        outfile = os.path.normpath(outfile)
        outfile = os.path.realpath(outfile)
        print(outfile)

    # 2) Construct Stream object
    st = read(infile)

    # 3) Plot
    if outdayplot:
        plot_option_type = 'dayplot'
    else:
        plot_option_type = 'normal'
    if outfile is not None:
        outfile = str(outfile)
        # outfile = outfile.replace(".MSEED", "png")
        # outfile += '.png'
        if outfilter is not None:
            st.filter("lowpass", freq=int(outfilter), corners=10)   # , zerophase=True
        st.plot(type=plot_option_type, outfile=outfile, size=(800, 600))
    else:
        print(st[0].stats)
        print(st[0].data)
        if outfilter is not None:
            st.filter("lowpass", freq=int(outfilter), corners=10)   # , zerophase=True
        st.plot(type=plot_option_type, method='full')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file')
    parser.add_argument('--dayplot', action='store_true', help='dayplot of the given file(s), normally same channel')
    parser.add_argument('--filter', action='store', help='Filter signal before ploting')
    args = parser.parse_args()

    plot_file(infile=args.infile,
              outfile=args.outfile,
              outfilter=args.filter,
              outdayplot=args.dayplot)
