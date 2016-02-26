#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from obspy.core import Trace, Stats, Stream, UTCDateTime
import numpy as np

parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
parser.add_argument('--outfile', action='store', help='name of output file', required=True)
parser.add_argument('--type', action='store', help='Type of output file', required=True)
args = parser.parse_args()

# 1) Make sure user inputs are correct
outfile_name = args.outfile
type_name = args.type
# Convert to real (no symlink) and full path
infile_paths = args.infile
for i in range(0, len(infile_paths)):
    infile_paths[i] = os.path.normcase(infile_paths[i])
    infile_paths[i] = os.path.normpath(infile_paths[i])
    infile_paths[i] = os.path.realpath(infile_paths[i])
    print(infile_paths[i])

# 1.5) If TUPLE, convert to Stream and Traces
if ".slist" in infile_paths[0]:
    #print("Tuple file %s" % infile_paths[0])
    print("SLIST file %s" % infile_paths[0])
    print("Creating MSEED files for every channel ..")

    fd = open(infile_paths[0])
    sl_header = fd.readline()
    sl_header = sl_header.replace("\r", "")
    sl_header = sl_header.replace("\n", "")
    sl_header = sl_header.replace("  ", "")
    sl_header = sl_header.replace(" ", "")
    sl_header = sl_header.split(",")
    sl_header_start_time = sl_header[3]
    print(sl_header)
    print(sl_header_start_time)

    #samples = fd.read()
    samples = []
    for line in fd:
        line = line.replace("\r", "")
        line = line.replace("\n", "")
        samples.append(line)
    fd.close()

    print("samples")
    print(samples)
    print("samples[0]")
    print(samples[0])
    fd.close()

    #sl_data = samples
    data = np.int32(samples)

    # # Convert to NumPy character array
    # data = np.fromstring(sl_data, dtype='|S1')
    # print(data)

    # Fill header attributes
    stats = {'network': 'BW', 'station': 'RJOB', 'location': 'ASD',
             'channel': 'SHZ', 'npts': len(data), 'sampling_rate': 100,
             'mseed': {'dataquality': 'D'},
             'starttime': UTCDateTime(str(sl_header_start_time))}
    # set current time
    #stats['starttime'] = UTCDateTime()
    st = Stream([Trace(data=data, header=stats)])
    # write as ASCII file (encoding=0)
    #outfile_name += ".mseed"
    #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
    st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')

    # Show that it worked, convert NumPy character array back to string
    st1 = read(outfile_name)
    print(st1[0].data.tostring())

    exit(0)

# 2) Construct Stream object
st = read(infile_paths[0])
for i in range(1, len(infile_paths)):
    st += read(infile_paths[i])
st = st.sort()

st[0].stats['mseed'] = {'dataquality': 'D'}
st[0].data = np.int32(st[0].data)

# 3) Convert
if type_name == "SLIST" or type_name == "slist":
    outfile_name += ".slist"
    st.write(outfile_name, format="SLIST")
    exit(0)
if type_name == "TSPAIR" or type_name == "tspair":
    outfile_name += ".tspair"
    st.write(outfile_name, format="TSPAIR")
    exit(0)
if type_name == "MSEED" or type_name == "mseed":
    outfile_name += ".mseed"
    st.write(outfile_name, format="MSEED")
    exit(0)