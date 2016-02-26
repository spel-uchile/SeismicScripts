#!/usr/bin/env python
__author__ = 'toopazo'

# Following tutorial at http://docs.obspy.org/tutorial/index.html

# 2) UTCDateTime
# from obspy.core import UTCDateTime
# time1 = UTCDateTime("2012-09-07T12:15:00")
# print(time1)
# print(time1.year)
# print(time1.julday)
# print(time1.timestamp)
# print(time1.weekday)
#
# time2 = UTCDateTime(2012, 1, 1)
# print(time2)
# print(time1 - time2)
# print(UTCDateTime())

# 3) Reading Seismograms
# from obspy.core import read
# #st = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')
# #st = read('seismonux/DATALOG/2015-09-23-1341-16S.R0023_001_SHZ')
# st1 = read('SPQRUSB/stream/2015-09-23T14-00-00-000Z_RST0000_0_ch4.MSEED')
# print(st1)
# print(len(st1))
# tr = st1[0]  # assign first and only trace to new variable
# print(tr)
# print(tr.stats)
# print(len(tr.data))
# print(tr.data)
# st1.plot()

# 4) Waveform Plotting Tutorial
# from obspy.core import read
# #singlechannel = read('http://examples.obspy.org/COP.BHZ.DK.2009.050')
# singlechannel = read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch4.MSEED')
# print(singlechannel)
# singlechannel.plot(outfile='singlechannel.png')
# #singlechannel.plot(type='dayplot')
#
# # threechannels = read('http://examples.obspy.org/COP.BHE.DK.2009.050')
# # threechannels += read('http://examples.obspy.org/COP.BHN.DK.2009.050')
# # threechannels += read('http://examples.obspy.org/COP.BHZ.DK.2009.050')
# threechannels = read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch4.MSEED')
# threechannels += read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch3.MSEED')
# threechannels += read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch2.MSEED')
# print(threechannels)
# threechannels.plot(outfile='threechannels.png')
# threechannels.plot(size=(800, 600))

# 6) Filtering Seismograms
# from obspy.core import read
# #st = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')
# #st = read('seismonux/DATALOG/2015-09-23-1341-16S.R0023_001_SHZ')
# st1 = read('SPQRUSB/stream/2015-09-23T14-00-00-000Z_RST0000_0_ch4.MSEED')
# print(st1)
# print(len(st1))
# tr = st1[0]  # assign first and only trace to new variable
# print(tr)
# print(tr.stats)
# print(len(tr.data))
# print(tr.data)
# st2 = st1.copy()
# st2.filter("lowpass", freq=7, corners=10)   # , zerophase=True
# st3 = st1
# st3 += st2
# print(st3)
# st3.plot()
# st3.plot(outfile='2015-09-23T14-00-00-000Z_RST0000_0_ch4_filtered.png')

# # 8) Merging Seismograms
# from obspy.core import read
# #st = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')
# #st = read('seismonux/DATALOG/2015-09-23-1341-16S.R0023_001_SHZ')
# st1 = read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch4.MSEED')
# tr1 = st1[0]  # assign first and only trace to new variable
# print(st1)
# #print(st1[0].stats)
# st1.plot(method='full')
# print("******************************")
#
# st2 = read('SPQRUSB/stream/2015-09-24T15-00-00-000Z_RST0000_0_ch4.MSEED')
# tr2 = st2[0]  # assign first and only trace to new variable
# print(st2)
# #print(st2[0].stats)
# st2.plot(method='full')
# print("******************************")
#
# st3 = st1
# st3 += st2
# print(st3)
# #print(st3[0].stats)
# #print(st3[1].stats)
# st3.plot(method='full')
# print("******************************")
#
# st3.merge()
# print(st3)
# #print(st3[0].stats)
# #print(st3[1].stats)
# st3.plot(method='full')
# st3.plot(outfile='2015-09-24T14-00-00-000Z_RST0000_0_ch4_merged.png')
# print("******************************")
#
# # Slicing Seismigrams
# from obspy.core import UTCDateTime
# dt1 = UTCDateTime("2015-09-24T14:00:00")
# dt2 = UTCDateTime("2015-09-24T15:00:00")
# st1 = st3.slice(dt1, dt2)
# st2 = st3.slice(dt2)
# print(st3)
# #print(st3[0].stats)
# #print(st3[1].stats)
# st1.plot(method='full')
# st2.plot(method='full')
# print("******************************")

# 11) Plotting Spectrograms
# from obspy.core import read
# #st = read("http://examples.obspy.org/RJOB_061005_072159.ehz.new")
# st = read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch4.MSEED')
# print(st)
# print(st[0])
# tr0_header = str(st[0])[:15]
# st.spectrogram(log=True, title=tr0_header + '    ' + str(st[0].stats.starttime),
#                outfile='2015-09-24T14-00-00-000Z_RST0000_0_ch4_spectrum.png')

# 12) Trigger/Picker Tutorial
from obspy.core import read
#from obspy.signal.trigger import plotTrigger
#from obspy.signal.trigger import classicSTALTA
#from obspy.signal.trigger import pkBaer
from obspy.signal.trigger import arPick
from obspy.core import UTCDateTime

#st = read("http://examples.obspy.org/ev0_6.a01.gse2")
st = read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch4.MSEED')
st += read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch3.MSEED')
st += read('SPQRUSB/stream/2015-09-24T14-00-00-000Z_RST0000_0_ch2.MSEED')
print(st)

tr1 = st[0]
#print(tr1.stats)
df = tr1.stats.sampling_rate
start_time = tr1.stats['starttime']
#print(start_time)
start_time = UTCDateTime(start_time)
#print(start_time)
#st.filter("lowpass", freq=7, corners=10)   # , zerophase=True
tr1 = st[0]
tr2 = st[1]
tr3 = st[2]

#p_pick, phase_info = pkBaer(tr1.data, df, 20, 60, 7.0, 12.0, 100, 100)
#print(p_pick)
#print(phase_info)
#print(p_pick/phase_info)
#st = st.slice(starttime=start_time, endtime=(start_time+15*60))

p_pick, s_pick = arPick(tr1.data, tr2.data, tr3.data, df, 1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)
print(p_pick)
print(s_pick)
p_pick = start_time + p_pick
s_pick = start_time + s_pick
print(p_pick)
print(s_pick)

st.plot()