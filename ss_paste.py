#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
# from obspy.core import UTCDateTime


def concat(infile, outfile):
    # 1) Make sure user inputs are correct
    # Convert to real (no symlink) and full path
    for i in range(0, len(infile)):
        infile[i] = os.path.normcase(infile[i])
        infile[i] = os.path.normpath(infile[i])
        infile[i] = os.path.realpath(infile[i])
        print(infile[i])
    # otufile
    outfile = str(outfile)
    outfile = os.path.normcase(outfile)
    outfile = os.path.normpath(outfile)
    outfile = os.path.realpath(outfile)

    # 2) Construct Stream object
    st = read(infile[0])
    for i in range(1, len(infile)):
        st += read(infile[i])
    st = st.sort()

    # 3) Save Stream
    # st = st.slice(starttime=start_time, endtime=stop_time)
    #st.write(filename=outfile, format='MSEED', byteorder='>')
    st.write(filename=outfile, format='MSEED')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Slices/cuts given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
    #parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('--start_time', action='store', help='start time (ISO 8601 format)', required=True)
    # parser.add_argument('--stop_time', action='store', help='stop time (ISO 8601 format)', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file', required=True)
    args = parser.parse_args()

    concat(infile=args.infile, outfile=args.outfile)