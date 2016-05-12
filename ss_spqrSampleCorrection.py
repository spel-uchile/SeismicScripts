#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
# import movingaverage


def samplecorrection(infile, outfile):
    # 1) Make sure user inputs are correct
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

    # # 3) Correct samples
    # for iter in range(0, 5):
    #     # samples = st[0].detrend().data
    #     averages = list(movingaverage.movingaverage(samples, 2, data_is_list=None, avoid_fp_drift=False))
    #     # arg = "len(samples) = %s =>len(averages) = %s " % (len(samples), len(averages))
    #     # print(arg)
    #     aver_prev = averages[0]
    #     for i in range(0, 182):
    #         arg = "samples[%s] = %s | averages[%s] = %s | aver_prev = %s" \
    #               % (i, samples[i], i, averages[i], aver_prev)
    #         print(arg)
    #         aver_prev = averages[i]
    #     aver_prev = averages[0]
    #     delta = 500
    #     for i in range(0, len(samples)-1):
    #         if not (abs(aver_prev) - delta) <= abs(averages[i]) <= (abs(aver_prev) + delta):
    #             try:
    #                 if samples[i] == -1804:
    #                     arg = "[correction] samples[%s] = %s and samples[%s] = %s | averages[%s] = %s aver_prev = %s"\
    #                           % (i, samples[i], i+1, samples[i+1], i, averages[i], aver_prev)
    #                     print(arg)
    #                 samples[i+1] = int(aver_prev)
    #                 averages[i] = aver_prev
    #                 averages[i+1] = aver_prev
    #                 if samples[i] == -1804:
    #                     arg = "[correction] samples[%s] = %s and samples[%s] = %s | averages[%s] = %s aver_prev = %s"\
    #                           % (i, samples[i], i+1, samples[i+1], i, averages[i], aver_prev)
    #                     print(arg)
    #             except IndexError:
    #                 pass
    #         else:
    #             aver_prev = averages[i]
    #     print(iter)

    # 3) Correct samples
    arg = "[samplecorrection] Correcting samples for %s .." % infile
    print(arg)
    samples = st[0].data
    delta = 1000
    for i in range(0, len(samples)):
        try:
            width = 2
            mean = 0
            for j in range(1, width+1):
                mean += samples[i-j]
            mean /= width
            # mean = samples[i-1]
            if not (mean - delta) <= samples[i] <= (mean + delta):
                # arg = "[correction] mean = %s | samples[%s] = %s | samples[%s] = %s | samples[%s] = %s"\
                #       % (mean, i-2, samples[i-2], i-1, samples[i-1], i, samples[i])
                # arg = "[correction] mean = %s | samples[%s] = %s" % (mean, i, samples[i])
                # print(arg)
                samples[i] = mean
        except IndexError:
            pass

    # 3) Save samples to MSEED
    arg = "[samplecorrection] Saving corrected file %s .." % infile
    print(arg)
    st[0].data = samples
    #st[0].data = st[0].detrend().data
    st[0].filter("lowpass", freq=20, corners=10)
    st.write(filename=outfile, format='MSEED')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file', required=True)
    args = parser.parse_args()

    samplecorrection(infile=args.infile, outfile=args.outfile)