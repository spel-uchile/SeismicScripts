#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np
import dateutil.parser


def separate_channels(infile):

    st = read(infile)
    len_st = len(st)
    if len_st == 1:
        arg = "[separate_channels] There is only one Trace in file file %s, nothing to do here .." % infile
        print(arg)
        return
    else:
        arg = "[separate_channels] There are %s Traces in file file %s, processing .." % (len_st, infile)
        print(arg)
        for tr_i in st:
            print(tr_i.stats)
            outfile_channel_ext = "_" + tr_i.stats.channel + ".MSEED"
            outfile = infile.replace(".MSEED", outfile_channel_ext)
            arg = "[separate_channels] Writing file %s .." % outfile
            print(arg)
            stc = tr_i
            stc.write(outfile, format='MSEED')
            # st.write(outfile, format='MSEED', encoding=11, reclen=256, byteorder='>')
            #print(tr_i.stats)
            print("")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', required=True)
    args = parser.parse_args()

    separate_channels(infile=args.infile)
