#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from obspy.core import UTCDateTime


def cut_stream(infile, starttime, endtime, outfile):
    # 1) Make sure user inputs are correct
    # Convert to real (no symlink) and full path
    infile = os.path.normcase(infile)
    infile = os.path.normpath(infile)
    infile = os.path.realpath(infile)
    print(infile)
    # Convert start and stop times to UTCDateTime
    starttime = UTCDateTime(starttime)
    print(starttime)
    endtime = UTCDateTime(endtime)
    print(endtime)
    # otufile
    if outfile is not None:
        outfile = os.path.normcase(outfile)
        outfile = os.path.normpath(outfile)
        outfile = os.path.realpath(outfile)
    else:
        outfile = str(starttime)
        outfile = outfile.replace(':', '-')
        outfile += '.MSEED'
    print(outfile)

    # 2) Construct Stream object
    st = read(infile)

    # 3) Slice Stream and save it
    st = st.slice(starttime=starttime, endtime=endtime)
    #st.write(filename=outfile_path, format='MSEED', byteorder='>')
    st.write(filename=outfile, format='MSEED')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Slices/cuts given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
    #parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    parser.add_argument('--starttime', action='store', help='start time (ISO 8601 format)', required=True)
    parser.add_argument('--endtime', action='store', help='stop time (ISO 8601 format)', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file')    # , required=True)
    args = parser.parse_args()

    cut_stream(infile=args.infile,
               starttime=args.starttime,
               endtime=args.endtime,
               outfile=args.outfile)