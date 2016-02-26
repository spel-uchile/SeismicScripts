#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from obspy.core import UTCDateTime

parser = argparse.ArgumentParser(description='Slices/cuts given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--start_time', action='store', help='start time (ISO 8601 format)', required=True)
parser.add_argument('--stop_time', action='store', help='stop time (ISO 8601 format)', required=True)
parser.add_argument('--outfile', action='store', help='name of output file', required=True)
args = parser.parse_args()

# 1) Make sure user inputs are correct
# Convert to real (no symlink) and full path
infile_paths = args.infile
for i in range(0, len(infile_paths)):
    infile_paths[i] = os.path.normcase(infile_paths[i])
    infile_paths[i] = os.path.normpath(infile_paths[i])
    infile_paths[i] = os.path.realpath(infile_paths[i])
    print(infile_paths[i])
# Convert start and stop times to UTCDateTime
start_time = str(args.start_time)
start_time = UTCDateTime(start_time)
print(str(args.start_time), start_time)
stop_time = str(args.stop_time)
stop_time = UTCDateTime(stop_time)
print(str(args.stop_time), stop_time)
# otufile
outfile_path = str(args.outfile)

# 2) Construct Stream object
st = read(infile_paths[0])
for i in range(1, len(infile_paths)):
    st += read(infile_paths[i])
st = st.sort()

# 3) Slice Stream and save it
st = st.slice(starttime=start_time, endtime=stop_time)
#st.write(filename=outfile_path, format='MSEED', byteorder='>')
st.write(filename=outfile_path, format='MSEED')