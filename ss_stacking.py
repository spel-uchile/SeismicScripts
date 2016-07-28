#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
import numpy as np


def stacking(infile1, infile2, outfile):
    # 1) Make sure user inputs are correct
    # Convert to real (no symlink) and full path
    infile1 = os.path.normcase(infile1)
    infile1 = os.path.normpath(infile1)
    infile1 = os.path.realpath(infile1)
    print(infile1)
    infile2 = os.path.normcase(infile2)
    infile2 = os.path.normpath(infile2)
    infile2 = os.path.realpath(infile2)
    print(infile1)
    if outfile is not None:
        outfile = os.path.normcase(outfile)
        outfile = os.path.normpath(outfile)
        outfile = os.path.realpath(outfile)
    else:
        # ss_utilities.ParserUtilities.
        infile1 = infile1.split("/")
        infile1 = infile1[len(infile1)-1]
        print(infile1)
        infile2 = infile2.split("/")
        infile2 = infile2[len(infile2)-1]
        print(infile2)
        outfile = infile1.replace(".MSEED", "_and_") + infile2.replace(".MSEED", "_stacking.MSEED")
    print(outfile)

    # 2) Construct Stream object
    st1 = read(infile1)
    stats1 = st1[0].stats
    st2 = read(infile2)
    stats2 = st2[0].stats

    # 3) Match starting and endtimes
    arg = "stats1 %s " % stats1
    print(arg)
    starttime1 = stats1['starttime']
    arg = "starttime1 %s " % starttime1
    print(arg)
    endtime1 = stats1['endtime']
    arg = "endtime1 %s " % endtime1
    print(arg)
    sampling_rate1 = stats1['sampling_rate']
    arg = "sampling_rate1 %s " % sampling_rate1
    print(arg)

    arg = "stats2 %s " % stats2
    print(arg)
    starttime2 = stats2['starttime']
    arg = "starttime2 %s " % starttime2
    print(arg)
    endtime2 = stats2['endtime']
    arg = "endtime2 %s " % endtime2
    print(arg)
    sampling_rate2 = stats2['sampling_rate']
    arg = "sampling_rate2 %s " % sampling_rate2
    print(arg)

    if sampling_rate1 != sampling_rate2:
        arg = "sampling_rate1 != sampling_rate2"
        print(arg)
        exit(0)
    if starttime1 >= endtime2:
        arg = "starttime1 >= endtime2"
        print(arg)
        exit(0)
    if starttime2 >= endtime1:
        arg = "starttime2 >= endtime1"
        print(arg)
        exit(0)

    if starttime1 > starttime2:
        starttime3 = starttime1
    else:
        starttime3 = starttime2
    if endtime1 > endtime2:
        endtime3 = endtime2
    else:
        endtime3 = endtime1
    arg = "starttime3 %s" % starttime3
    print(arg)
    arg = "endtime3 %s" % endtime3
    print(arg)
    st1 = st1.slice(starttime=starttime3, endtime=endtime3)
    samples1 = np.float64(st1[0].data)
    st2 = st2.slice(starttime=starttime3, endtime=endtime3)
    samples2 = np.float64(st2[0].data)

    # 3) Stack commmon signals
    samples3 = (samples1 + samples2)/2
    # arg = ""
    # print(arg)

    # Write to MSEED
    st1[0].data = np.int32(samples3)
    st1.write(outfile, format='MSEED', encoding=11, reclen=256, byteorder='>')
    #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
    arg = "[tuple2mseed] MSEED created: %s" % outfile
    print(arg)
    # st1 = read(outfile)
    # arg = "[tuple2mseed] MSEED created: %s" % st1[0]
    # print(arg)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile1', action='store', help='files to process', required=True)
    parser.add_argument('--infile2', action='store', help='files to process', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file')
    args = parser.parse_args()

    stacking(infile1=args.infile1,
             infile2=args.infile2,
             outfile=args.outfile)
