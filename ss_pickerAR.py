#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
#from obspy.core import UTCDateTime
from obspy.signal.trigger import arPick
#from obspy.signal.trigger import pkBaer
from obspy.core import UTCDateTime

parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
parser.add_argument('--infile', action='store', help='files to process', nargs='+', required=True)
#parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--outfile', action='store', help='name of output file')    # , required=True)
# parser.add_argument('--dayplot', action='store_true', help='dayplot of the given file(s), normally same channel')
args = parser.parse_args()

# 1) Make sure user inputs are correct
outfile_plot = args.outfile
# dayplot = args.dayplot
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

# 3) Do the picking
tr1 = st[0]
#print(tr1.stats)
df = tr1.stats.sampling_rate
start_time = tr1.stats['starttime']
# #print(start_time)
# start_time = UTCDateTime(start_time)
# #print(start_time)
# #st.filter("lowpass", freq=7, corners=10)   # , zerophase=True
tr1 = st[0]
tr2 = st[1]
tr3 = st[2]

# p_pick, phase_info = pkBaer(tr1.data, df, 20, 60, 7.0, 12.0, 100, 100)
# print(p_pick)
# print(phase_info)
# print(p_pick/df)
# print(start_time + (p_pick/df))
#st = st.slice(starttime=start_time, endtime=(start_time+15*60))

p_pick, s_pick = arPick(tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
print(p_pick)
print(s_pick)
p_pick = start_time + p_pick
s_pick = start_time + s_pick
print(p_pick)
print(s_pick)

#st.plot(type='relative')
#st.plot()
# outfile_plot += ".txt"
# fd = open(outfile_plot, "w")
# line = str(start_time + (p_pick/df))
# fd.write(line)
# fd.close()