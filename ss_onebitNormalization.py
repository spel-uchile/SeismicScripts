#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
# from obspy.core import Trace, Stream, UTCDateTime
# import numpy as np
# import dateutil.parser


def onebit_normalization(infile,  outfile):
    # 1) Make sure user inputs are correct (Convert to real -no symlink- and full path)
    infile = os.path.normcase(infile)
    infile = os.path.normpath(infile)
    infile = os.path.realpath(infile)
    print(infile)
    outfile = os.path.normcase(outfile)
    outfile = os.path.normpath(outfile)
    outfile = os.path.realpath(outfile)
    print(outfile)

    # 2) Open file and get data
    st1 = read(infile)
    arg = "[onebit_normalization] MSEED created: %s" % st1[0]
    print(arg)
    print(st1[0])
    # print(st1[0].stats)
    # print(st1[0].data)

    #st1[0].data = st1[0].detrend().data
    samples = st1[0].data
    mean = st1[0].data.mean()
    arg = "[onebit_normalization] mean %s" % mean
    print(arg)

    for i in range(0, len(samples)):
        if samples[i] >= mean:
            samples[i] = +1
        else:
            samples[i] = -1
    st1[0].data = samples

    mean = st1[0].data.mean()
    arg = "[onebit_normalization] mean %s" % mean
    print(arg)
    # st1[0].detrend()

    # ) Write data
    arg = "[onebit_normalization] writing %s .." % outfile
    print(arg)
    st1.write(outfile, format='MSEED', encoding=11, reclen=256, byteorder='>')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', required=True)
    parser.add_argument('--outfile', action='store', help='files to process', required=True)
    args = parser.parse_args()

    onebit_normalization(infile=args.infile, outfile=args.outfile)
