#!/usr/bin/env python
__author__ = 'toopazo'

import argparse
import os
from obspy.core import read
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np
import dateutil.parser


class TupleParser():
    """ Handler Class that takes a ".tuple" file and parse it. Detect errors, correct them, etc """
    def __init__(self):
        pass

    @staticmethod
    def get_data_and_correct_missing_samples(infile_path, starttime, sps):
        tuple_path = infile_path
        filename_date = starttime
        header_sps = sps

        matrix = TupleParser.load_file_to_matrix(pathfile=tuple_path, ncol=5, skipnrows=2)
        # matrix = TupleParser.discard_oversamples(matrix)

        # 4) Check correct SPS over all the file content and re-pack into Tuple object
        timestamps = []
        ch1_samples = []
        ch2_samples = []
        ch3_samples = []
        ch4_samples = []

        # New format
        # 15:00:00.010000 0 -2003 1371 -46306
        #datetime.strptime(matrix[0][0], "%H:%M:%S.%fZ")
        current_sec = dateutil.parser.parse(matrix[0][0])
        # Append samples to channels
        val_datetime = matrix[0][0]
        timestamps.append(val_datetime)
        ch1_samples.append(matrix[0][1])
        ch2_samples.append(matrix[0][2])
        ch3_samples.append(matrix[0][3])
        ch4_samples.append(matrix[0][4])

        sample_cnt = 1
        issues_cnt = 0
        print "[TupleParser.get_data_and_correct_missing_samples] current_sec %s " % current_sec
        print "[TupleParser.get_data_and_correct_missing_samples] sample_cnt %s " % sample_cnt
        for i in range(1, len(matrix)):
            sec_i = dateutil.parser.parse(matrix[i][0])
            sample_cnt += 1
            # print "sec_i %s " % sec_i
            if sec_i.second != current_sec.second:
                sample_cnt -= 1     # discount extra sample (from next second)
                if sample_cnt == header_sps:
                    # print "[TupleParser.get_data_and_correct_missing_samples] ***************************************"
                    # print "[TupleParser.get_data_and_correct_missing_samples] %s != %s [sec_i | current_sec] " %\
                    #       (sec_i.second, current_sec.second)
                    # print "[TupleParser.get_data_and_correct_missing_samples] %s == %s [sample_cnt | header_sps]" %\
                    #       (sample_cnt, self.header_sps)
                    # print "[TupleParser.get_data_and_correct_missing_samples] sec_i %s " % sec_i
                    # print "[TupleParser.get_data_and_correct_missing_samples] current_sec %s " % current_sec
                    # print "[TupleParser.get_data_and_correct_missing_samples] sample_cnt %s " % sample_cnt
                    # print "[TupleParser.get_data_and_correct_missing_samples] matrix[%s] %s " % (i, matrix[i])
                    pass
                elif sample_cnt == header_sps - 1:
                    # print "[TupleParser.get_data_and_correct_missing_samples] ***************************************"
                    # print "[TupleParser.get_data_and_correct_missing_samples] %s != %s [sec_i | current_sec] " % \
                    #       (sec_i.second, current_sec.second)
                    # print "[TupleParser.get_data_and_correct_missing_samples] %s == %s - 1 [sample_cnt|header_sps-1]"\
                    #       % (sample_cnt, header_sps)
                    # print "[TupleParser.get_data_and_correct_missing_samples] sec_i %s " % sec_i
                    # print "[TupleParser.get_data_and_correct_missing_samples] current_sec %s " % current_sec
                    # print "[TupleParser.get_data_and_correct_missing_samples] sample_cnt %s " % sample_cnt
                    # print "[TupleParser.get_data_and_correct_missing_samples] matrix[%s] %s " % (i, matrix[i])

                    # Add sample number 100 (repeat previous sample)
                    # print "[TupleParser.get_data_and_correct_missing_samples] Adding missing sample .."
                    last_sec = dateutil.parser.parse(matrix[i-1][0])
                    extra_time = last_sec + (sec_i - last_sec)/2
                    extra_time = str(extra_time)[11:]
                    val_datetime = "%sT%s" % (filename_date, extra_time)
                    # print "[TupleParser.get_data_and_correct_missing_samples] val_datetime %s " % val_datetime
                    timestamps.append(val_datetime)
                    val = (int(matrix[i][1]) + int(matrix[i-1][1]))/2
                    # print "[TupleParser.get_data_and_correct_missing_samples] ch1 %s " % val
                    ch1_samples.append(val)
                    val = (int(matrix[i][2]) + int(matrix[i-1][2]))/2
                    # print "[TupleParser.get_data_and_correct_missing_samples] ch2 %s " % val
                    ch2_samples.append(val)
                    val = (int(matrix[i][3]) + int(matrix[i-1][3]))/2
                    # print "[TupleParser.get_data_and_correct_missing_samples] ch3 %s " % val
                    ch3_samples.append(val)
                    val = (int(matrix[i][4]) + int(matrix[i-1][4]))/2
                    # print "[TupleParser.get_data_and_correct_missing_samples] ch4 %s " % val
                    ch4_samples.append(val)

                    issues_cnt += 1
                else:
                    print "[TupleParser.get_data_and_correct_missing_samples] *****************************************"
                    print "[TupleParser.get_data_and_correct_missing_samples] %s != %s [sec_i | current_sec] " %\
                          (sec_i.second, current_sec.second)
                    print "[TupleParser.get_data_and_correct_missing_samples] %s != %s - 1 [sample_cnt | header_sps-1]"\
                          % (sample_cnt, header_sps)
                    print "[TupleParser.get_data_and_correct_missing_samples] sec_i %s " % sec_i
                    print "[TupleParser.get_data_and_correct_missing_samples] current_sec %s " % current_sec
                    print "[TupleParser.get_data_and_correct_missing_samples] sample_cnt %s " % sample_cnt
                    print "[TupleParser.get_data_and_correct_missing_samples] matrix[%s] %s " % (i, matrix[i])
                    issues_cnt += 1

                # Update counters
                current_sec = sec_i
                sample_cnt = 1

            # Append samples to channels
            val_datetime = "%sT%s" % (filename_date, matrix[i][0])
            timestamps.append(val_datetime)
            ch1_samples.append(matrix[i][1])
            ch2_samples.append(matrix[i][2])
            ch3_samples.append(matrix[i][3])
            ch4_samples.append(matrix[i][4])

        print "[TupleParser.get_data_and_correct_missing_samples] *****************************************************"
        print "[TupleParser.get_data_and_correct_missing_samples] issues_cnt %s (missing samples) " % issues_cnt

        # display a few samples
        print "[TupleParser.get_data_and_correct_missing_samples] matrix[0] = %s" % matrix[0]
        print "[TupleParser.get_data_and_correct_missing_samples] matrix[1] = %s" % matrix[1]
        print "[TupleParser.get_data_and_correct_missing_samples] matrix[2] = %s" % matrix[2]
        print "[TupleParser.get_data_and_correct_missing_samples] *****************************************************"

        resp_dict = {'time': timestamps,
                     'ch4': ch4_samples,
                     'ch3': ch3_samples,
                     'ch2': ch2_samples,
                     'ch1': ch1_samples}
        return resp_dict

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
        resp_dict = TupleParser.get_data_and_correct_missing_samples(infile_path=infile,
                                                                     starttime=tuple_header_starttime,
                                                                     sps=tuple_header_sps)
        # time = resp_dict['time']
        data_ch4 = resp_dict['ch4']
        data_ch3 = resp_dict['ch3']
        data_ch2 = resp_dict['ch2']
        data_ch1 = resp_dict['ch1']

        # 5) Convert channel ch1
        data = np.int32(data_ch1)
        channel = user_ch1
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])
        outfile_name = tuple_header_starttime + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

        # 6) Convert channel SHZ
        data = np.int32(data_ch2)
        channel = user_ch2
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])
        outfile_name = tuple_header_starttime + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

        # 7) Convert channel SHE
        data = np.int32(data_ch3)
        channel = user_ch3
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])
        outfile_name = tuple_header_starttime + "_" + tuple_filename_station + "_" + channel + ".MSEED"
        outfile_name = outfile_name.replace(":", "-")
        st.write(outfile_name, format='MSEED', encoding=11, reclen=256, byteorder='>')
        #st.write(outfile_name, format='MSEED', encoding=0, reclen=256)
        st1 = read(outfile_name)
        arg = "[tuple2mseed] MSEED created: %s" % st1[0]
        print(arg)
        # print(st1[0])
        # print(st1[0].stats)
        # print(st1[0].data)

        # 8) Convert channel ch4
        data = np.int32(data_ch4)
        channel = user_ch4
        # Fill header attributes
        stats = {'network': 'UNK', 'station': tuple_filename_station, 'location': 'UNK',
                 'channel': channel, 'npts': len(data), 'sampling_rate': tuple_header_sps,
                 'mseed': {'dataquality': 'D'},
                 'starttime': UTCDateTime(str(tuple_header_starttime))}
        st = Stream([Trace(data=data, header=stats)])
        outfile_name = tuple_header_starttime + "_" + tuple_filename_station + "_" + channel + ".MSEED"
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
