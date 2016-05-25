#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
# from obspy.core import read, Stream, Trace
import numpy as np


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
    else:
        outfile = infile
    print(outfile)

    # 2) Construct Stream object and iterate over every Trace in Stream
    st = read(infile)
    for tr_i in range(0, len(st)):

        # 3) Correct samples,  correct very high freq noise (signal "jumps")
        arg = "[samplecorrection] Correcting samples for %s .." % infile
        print(arg)
        samples = st[tr_i].data
        delta = 1000
        for sample_i in range(0, len(samples)):
            try:
                width = 2
                mean = 0
                for sample_j in range(1, width+1):
                    mean += samples[sample_i-sample_j]
                mean /= width
                # mean = samples[sample_i-1]
                if not (mean - delta) <= samples[sample_i] <= (mean + delta):
                    # arg = "[correction] mean = %s | samples[%s] = %s | samples[%s] = %s | samples[%s] = %s"\
                    #       % (mean, sample_i-2, samples[sample_i-2], sample_i-1,
                    #  samples[sample_i-1], sample_i, samples[sample_i])
                    # arg = "[correction] mean = %s | samples[%s] = %s" % (mean, sample_i, samples[sample_i])
                    # print(arg)
                    samples[sample_i] = mean
            except IndexError:
                pass
        st[tr_i].data = samples

        # 4) Substract average (self-made DC detrend, actually works better than obspy.detrend )
        samples = st[tr_i].data
        samples = np.int32(samples)
        samples_average = samples.mean()
        for sample_j in range(0, len(samples)):
            samples[sample_j] -= samples_average
        st[tr_i].data = samples
        # print(st[tr_i])
        # print(st[tr_i].data)

        # 5) Obspy filter and detrend
        #st[0].data = st[0].detrend().data
        st[tr_i].filter("lowpass", freq=40, corners=10)

        # End of loop activities

    arg = "[samplecorrection] Saving corrected file %s .." % infile
    print(arg)
    st.write(filename=outfile, format='MSEED')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', required=True)
    parser.add_argument('--outfile', action='store', help='name of output file', required=True)
    args = parser.parse_args()

    samplecorrection(infile=args.infile, outfile=args.outfile)


        # # 3) Correct samples
    # for iter in range(0, 5):
    #     # samples = st[0].detrend().data
    #     averages = list(movingaverage.movingaverage(samples, 2, data_is_list=None, avoid_fp_drift=False))
    #     # arg = "len(samples) = %s =>len(averages) = %s " % (len(samples), len(averages))
    #     # print(arg)
    #     aver_prev = averages[0]
    #     for sample_i in range(0, 182):
    #         arg = "samples[%s] = %s | averages[%s] = %s | aver_prev = %s" \
    #               % (sample_i, samples[sample_i], sample_i, averages[sample_i], aver_prev)
    #         print(arg)
    #         aver_prev = averages[sample_i]
    #     aver_prev = averages[0]
    #     delta = 500
    #     for sample_i in range(0, len(samples)-1):
    #         if not (abs(aver_prev) - delta) <= abs(averages[sample_i]) <= (abs(aver_prev) + delta):
    #             try:
    #                 if samples[sample_i] == -1804:
    #                     arg = "[correction] samples[%s] = %s and samples[%s] = %s | averages[%s] = %s aver_prev = %s"\
    #                           % (sample_i, samples[sample_i], sample_i+1,
    # samples[sample_i+1], sample_i, averages[sample_i], aver_prev)
    #                     print(arg)
    #                 samples[sample_i+1] = int(aver_prev)
    #                 averages[sample_i] = aver_prev
    #                 averages[sample_i+1] = aver_prev
    #                 if samples[sample_i] == -1804:
    #                     arg = "[correction] samples[%s] = %s and samples[%s] = %s | averages[%s] = %s aver_prev = %s"\
    #                           % (sample_i, samples[sample_i], sample_i+1,
    # samples[sample_i+1], sample_i, averages[sample_i], aver_prev)
    #                     print(arg)
    #             except IndexError:
    #                 pass
    #         else:
    #             aver_prev = averages[sample_i]
    #     print(iter)