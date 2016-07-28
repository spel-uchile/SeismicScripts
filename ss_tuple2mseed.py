#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np
import dateutil.parser
from datetime import datetime, timedelta


class TupleParser():
    """ Handler Class that takes a ".tuple" file and parse it. Detect errors, correct them, etc """
    def __init__(self):
        pass

    @staticmethod
    def get_and_correct_missing_samples(infile_path, starttime, sps):
        tuple_path = infile_path
        filename_date = starttime

        matrix = TupleParser.load_file_to_matrix(pathfile=tuple_path, ncol=5, skipnrows=2)
        # matrix = TupleParser.discard_oversamples(matrix)
        print "*****************************************************"
        print "[TupleParser.get_and_correct_missing_samples] Before correction:"
        print "[TupleParser.get_and_correct_missing_samples] matrix[0] = %s" % matrix[0]
        print "[TupleParser.get_and_correct_missing_samples] matrix[1] = %s" % matrix[1]
        print "[TupleParser.get_and_correct_missing_samples] matrix[2] = %s" % matrix[2]
        print "[TupleParser.get_and_correct_missing_samples] matrix[%s] = %s" % \
              (len(matrix) - 1, matrix[len(matrix) - 1])

        # 4) Check correct SPS over all the file content and re-pack into Tuple object
        timestamps = []
        ch1_samples = []
        ch2_samples = []
        ch3_samples = []
        ch4_samples = []

        print "*****************************************************"
        first_microsec_sps = timedelta(microseconds=0).microseconds
        print "first_microsec_sps %s" % first_microsec_sps
        second_microsec_normal_sps = timedelta(microseconds=(1000000/sps)).microseconds
        print "second_microsec_normal_sps %s" % second_microsec_normal_sps
        second_microsec_1missing_sps = timedelta(microseconds=(1000000/(sps-1))).microseconds
        print "second_microsec_1missing_sps %s" % second_microsec_1missing_sps
        second_microsec_1extra_sps = timedelta(microseconds=(1000000/(sps+1))).microseconds
        print "second_microsec_1extra_sps %s" % second_microsec_1extra_sps
        print "*****************************************************"
        sample_i = 0
        seconds_cnt = 0
        missing_sample_cnt = 0
        continuity_error_cnt = 0
        extra_sample_cnt = 0
        while True:
            # >> Check first sample
            first_timestamp = dateutil.parser.parse(matrix[sample_i][0])
            # print "Checking => first_timestamp %s" % first_timestamp

            # >> Check if all is right with the first sample
            if first_timestamp.microsecond != first_microsec_sps:
                print "Error at %s => %s != %s [first_timestamp.microsecond != first_microsec_sps]" %\
                      (first_timestamp, first_timestamp.microsecond, first_microsec_sps)
                raise RuntimeError

            # >> Append first sample
            val_datetime = "%sT%s" % (filename_date, matrix[sample_i][0])
            timestamps.append(val_datetime)
            ch1_samples.append(matrix[sample_i][1])
            ch2_samples.append(matrix[sample_i][2])
            ch3_samples.append(matrix[sample_i][3])
            ch4_samples.append(matrix[sample_i][4])

            # >> Check if second sample is of type "sps-1" or "sps" or "sps+1"
            sample_i += 1
            second_timestamp = dateutil.parser.parse(matrix[sample_i][0])
            # print "Checking => second_timestamp %s" % second_timestamp
            if second_timestamp.microsecond == second_microsec_normal_sps:
                # >> Append "second" sample
                val_datetime = "%sT%s" % (filename_date, matrix[sample_i][0])
                timestamps.append(val_datetime)
                ch1_samples.append(matrix[sample_i][1])
                ch2_samples.append(matrix[sample_i][2])
                ch3_samples.append(matrix[sample_i][3])
                ch4_samples.append(matrix[sample_i][4])
            elif second_timestamp.microsecond == second_microsec_1extra_sps:
                # >> Correct extra sample
                # .. do not append

                # >> Report issue
                extra_sample_cnt += 1
                # print "Correcting extra sample => extra_sample_cnt %s " % extra_sample_cnt
                # print "*****************************************************"
                pass
            elif second_timestamp.microsecond == second_microsec_1missing_sps:
                # >> Append missing sample as missing_sample = (second_sample + first_sample) / 2
                val = (int(matrix[sample_i][1]) + int(matrix[sample_i-1][1]))/2
                # print "[TupleParser.get_and_correct_missing_samples] ch1 %s " % val
                ch1_samples.append(val)
                val = (int(matrix[sample_i][2]) + int(matrix[sample_i-1][2]))/2
                # print "[TupleParser.get_and_correct_missing_samples] ch2 %s " % val
                ch2_samples.append(val)
                val = (int(matrix[sample_i][3]) + int(matrix[sample_i-1][3]))/2
                # print "[TupleParser.get_and_correct_missing_samples] ch3 %s " % val
                ch3_samples.append(val)
                val = (int(matrix[sample_i][4]) + int(matrix[sample_i-1][4]))/2
                # print "[TupleParser.get_and_correct_missing_samples] ch4 %s " % val
                ch4_samples.append(val)

                # >> Report issue
                missing_sample_cnt += 1
                # print "Correcting missing sample => missing_sample_cnt %s " % missing_sample_cnt
                # print "*****************************************************"

                # >> Append "second" sample
                val_datetime = "%sT%s" % (filename_date, matrix[sample_i][0])
                timestamps.append(val_datetime)
                ch1_samples.append(matrix[sample_i][1])
                ch2_samples.append(matrix[sample_i][2])
                ch3_samples.append(matrix[sample_i][3])
                ch4_samples.append(matrix[sample_i][4])
            else:
                print "Error at %s => %s, %s, %s (second_timestamp != sps) or (second_timestamp != sps-1)" \
                      " or (second_timestamp != sps+1) " %\
                      (second_timestamp, second_microsec_normal_sps,
                       second_microsec_1missing_sps, second_microsec_1extra_sps)
                raise RuntimeError

            # Append the rest of the samples in the current sec
            sample_i_timestamp = None
            while True:
                sample_i += 1
                try:
                    sample_i_timestamp = dateutil.parser.parse(matrix[sample_i][0])
                except IndexError:
                    resp_dict = {'ch4': ch4_samples,
                                 'ch3': ch3_samples,
                                 'ch2': ch2_samples,
                                 'ch1': ch1_samples}
                    print "len(matrix) = %s " % len(matrix)
                    print "missing_sample_cnt = %s " % missing_sample_cnt
                    print "continuity_error_cnt = %s " % continuity_error_cnt
                    print "seconds_cnt = %s " % seconds_cnt
                    print "sample_i = %s " % sample_i
                    print "len(ch1_samples) = %s " % len(ch1_samples)
                    print "len(ch2_samples) = %s " % len(ch2_samples)
                    print "len(ch3_samples) = %s " % len(ch3_samples)
                    print "len(ch4_samples) = %s " % len(ch4_samples)
                    print "*****************************************************"
                    return resp_dict
                if sample_i_timestamp.second != first_timestamp.second:
                    break
                # Append sample_i
                val_datetime = "%sT%s" % (filename_date, matrix[sample_i][0])
                timestamps.append(val_datetime)
                ch1_samples.append(matrix[sample_i][1])
                ch2_samples.append(matrix[sample_i][2])
                ch3_samples.append(matrix[sample_i][3])
                ch4_samples.append(matrix[sample_i][4])

            # >> Add current second to the counter
            seconds_cnt += 1

            # >> Check seconds continuity
            if first_timestamp + timedelta(seconds=1) != sample_i_timestamp:
                print "Error at %s => %s != %s (first_timestamp + timedelta(seconds=1) != sample_i_timestamp)" % \
                      (sample_i_timestamp, first_timestamp + timedelta(seconds=1), sample_i_timestamp)
                print "*****************************************************"
                continuity_error_cnt += 1
                # raise RuntimeError

    @staticmethod
    def get_info(infile_path):
        containing_folder, infile_name = os.path.split(infile_path)
        # arg = "head %s tail %s" % (containing_folder, infile_name)
        # print(arg)
        tuple_filename_station = infile_name.split("_")[1]
        arg = "[TupleParser.get_info] tuple_filename_station %s " % tuple_filename_station
        print(arg)

        fd = open(infile_path)
        fd_line1_header = fd.readline()
        # fd_line1_header = fd_line1_header.replace("\r", "")
        # fd_line1_header = fd_line1_header.replace("\n", "")
        # fd_line1_header = fd_line1_header.replace("  ", "")
        # fd_line1_header = fd_line1_header.replace(" ", "")
        # fd_line1_header = fd_line1_header.split(",")
        fd_line1_header = fd_line1_header.split(" ")
        tuple_header_sps = int(fd_line1_header[1])
        arg = "[TupleParser.get_info] tuple_header_sps %s" % tuple_header_sps
        print(arg)
        fd_line2_header = fd.readline()
        fd_line2_header.split(" ")
        fd_line3_header = fd.readline()
        fd_line3_header = fd_line3_header.split(" ")
        tuple_header_starttime = fd_line3_header[0]
        tuple_header_starttime = "%sT%s" % (infile_name[:10], tuple_header_starttime)
        arg = "[TupleParser.get_info] tuple_header_starttime %s" % tuple_header_starttime
        print(arg)

        #samples = fd.read()
        samples = []
        for line in fd:
            line = line.replace("\r", "")
            line = line.replace("\n", "")
            samples.append(line)
        fd.close()

        resp_dict = {'station': tuple_filename_station,
                     'sps': tuple_header_sps,
                     'starttime': tuple_header_starttime}
        return resp_dict

    @staticmethod
    def split_from_to(strg, s_from, s_to):
        arg = strg.split(s_from)
        #print(arg)
        arg = arg[1].split(s_to)
        #print(arg)
        arg = arg[0]
        #print(arg)
        return arg

    @staticmethod
    def load_file_to_matrix(pathfile, ncol, skipnrows):
        # Streamer_0, 100 SPS
        # GpsTime NotConnected SHE(EW) SHN(NS) SHZ(UD)
        # 02:00:00 1499 -2197 1844 -1881
        # 02:00:00.010000 1516 -2167 1849 -1758
        # 02:00:00.020000 1526 -2233 1834 -1815
        with open(pathfile, "r") as f:
            linecnt = 1
            lmatrix = []
            for f_line in f:
                if skipnrows >= linecnt:
                    pass
                else:
                    #clean EOL and preceeding characters
                    f_line = f_line.replace("      ", " ")
                    f_line = f_line.replace("     ", " ")
                    f_line = f_line.replace("    ", " ")
                    f_line = f_line.replace("   ", " ")
                    f_line = f_line.replace("  ", " ")
                    f_line = f_line.replace(" \r", "")
                    f_line = f_line.replace(" \n", "")
                    f_line = f_line.replace("\r", "")
                    f_line = f_line.replace("\n", "")

                    # split and calculate lenght
                    f_row = f_line.split(" ")
                    f_row_len = len(f_row)
                    if f_row[f_row_len-1] == "":
                        f_row = f_row[:-1]
                        f_row_len = len(f_row)

                    # save to matrix or report
                    if f_row_len != ncol:
                        print "[load_file_to_matrix] Error at line %s: ncol      : %s" % (linecnt, ncol)
                        print "[load_file_to_matrix] Error at line %s: f_line    : %s" % (linecnt, f_line)
                        print "[load_file_to_matrix] Error at line %s: f_row     : %s" % (linecnt, f_row)
                        print "[load_file_to_matrix] Error at line %s: f_row_len : %s" % (linecnt, f_row_len)
                        raise ValueError
                    else:
                        lmatrix.append(f_row)
                linecnt += 1
        return lmatrix

    @staticmethod
    def discard_oversamples(matrix):
        numrows = len(matrix)       # rows
        numcols = len(matrix[0])    # columns (samples)
        print("[discard_oversamples] numrows %s" % numrows)
        print("[discard_oversamples] numcols %s" % numcols)

        cu_sps = 0
        cu_sec = matrix[0][3]
        cu_sec = cu_sec[:8]
        for i in range(0, numrows):
            sec_i = matrix[i][3]
            sec_i = sec_i[:8]
            if sec_i == cu_sec:
                cu_sps += 1
                if cu_sps > 111:
                    print("[discard_oversamples] found cu_sps %s" % cu_sps)
                    matrix[i].pop(0)
                    matrix[i].pop(1)
                    matrix[i].pop(2)
                    matrix[i].pop(3)
                    matrix[i].pop(4)
                    matrix[i].pop(5)
            else:   # next second
                cu_sec = sec_i

        return matrix


def tuple2mseed(infile, user_ch1, user_ch2, user_ch3, user_ch4):
    # 1) Make sure user inputs are correct (Convert to real -no symlink- and full path)
    infile = os.path.normcase(infile)
    infile = os.path.normpath(infile)
    infile = os.path.realpath(infile)
    print(infile)

    # 2) If TUPLE, convert to MSEED
    if ".tuple" in infile:

        # 3) Get header and extra info
        arg = "[tuple2mseed] \".tuple\" file %s" % infile
        print(arg)
        resp_dict = TupleParser.get_info(infile)
        tuple_filename_station = resp_dict['station']
        tuple_header_sps = resp_dict['sps']
        tuple_header_starttime = resp_dict['starttime']
        # exit(0)

        # 4) Get samples
        arg = "[tuple2mseed] Creating MSEED files for every channel .."
        print(arg)
        resp_dict = TupleParser.get_and_correct_missing_samples(infile_path=infile,
                                                                starttime=tuple_header_starttime,
                                                                sps=tuple_header_sps)
        data_ch4 = resp_dict['ch4']
        data_ch3 = resp_dict['ch3']
        data_ch2 = resp_dict['ch2']
        data_ch1 = resp_dict['ch1']

        # >> Substract DC component
        data = np.int32(data_ch1)
        data = data - np.mean(data)
        data = np.int32(data)
        # >> Convert channel to MSEED
        channel = user_ch1
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])

        # >> Write to disk
        # print(tuple_header_starttime)
        outfile_name = tuple_header_starttime.split(".")
        # print(tuple_header_starttime)
        outfile_name = outfile_name[0]
        outfile_name = outfile_name + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

        # >> Substract DC component
        data = np.int32(data_ch2)
        data = data - np.mean(data)
        data = np.int32(data)
        # >> Convert channel to MSEED
        channel = user_ch2
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])

        # >> Write to disk
        # print(tuple_header_starttime)
        outfile_name = tuple_header_starttime.split(".")
        # print(tuple_header_starttime)
        outfile_name = outfile_name[0]
        outfile_name = outfile_name + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

        # >> Substract DC component
        data = np.int32(data_ch3)
        data = data - np.mean(data)
        data = np.int32(data)
        # >> Convert channel to MSEED
        channel = user_ch3
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])

        # >> Write to disk
        # print(tuple_header_starttime)
        outfile_name = tuple_header_starttime.split(".")
        # print(tuple_header_starttime)
        outfile_name = outfile_name[0]
        outfile_name = outfile_name + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

        # >> Substract DC component
        data = np.int32(data_ch4)
        data = data - np.mean(data)
        data = np.int32(data)
        # >> Convert channel to obspy Stream
        channel = user_ch4
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])

        # >> Write to disk in MSEED format
        # print(tuple_header_starttime)
        outfile_name = tuple_header_starttime.split(".")
        # print(tuple_header_starttime)
        outfile_name = outfile_name[0]
        outfile_name = outfile_name + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

    else:
        arg = "File %s does NOT end with \".tuple\"" % infile[0]
        print(arg)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='files to process', required=True)
    parser.add_argument('--ch1', action='store', help='channel signal', required=True)
    parser.add_argument('--ch2', action='store', help='channel signal', required=True)
    parser.add_argument('--ch3', action='store', help='channel signal', required=True)
    parser.add_argument('--ch4', action='store', help='channel signal', required=True)
    args = parser.parse_args()

    tuple2mseed(infile=args.infile,
                user_ch1=args.ch1,
                user_ch2=args.ch2,
                user_ch3=args.ch3,
                user_ch4=args.ch4)
