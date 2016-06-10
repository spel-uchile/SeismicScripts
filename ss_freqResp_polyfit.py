#!/usr/bin/env python

__author__ = 'toopazo'
 # # # -*- coding: utf-8 -*-

import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
import math
import cmath
import argparse
import ss_fft
import matplotlib.lines as mlines
import matplotlib.patches as mpatches


def evaluate_functions(eval_freq, polyfit_roots, polyfit_degrees, polyfit_coeffs):
    jw = 1j*(2.0*math.pi*eval_freq)

    # Polynomial fit
    # F(s = jw)  = w**n * p[0] + ... + w**1 * p[n-1] + p[n] = y(jw) [dB]
    yfit_val = 0.0
    for degree_j in range(0, len(polyfit_coeffs)):
        # x[0]**n * p[0] + ... + x[0] * p[n-1] + p[n] = y[0]
        yfit_val += (eval_freq**(polyfit_degrees-degree_j))*polyfit_coeffs[degree_j]
    # yfit_amp_val = abs(yfit_val)

    # GPZ Transfer Function
    # F(s = jw)  = (jw - z1) .. (jw - zn) / (jw - p1) .. (jw - pm)  = y(jw)
    ygpz_val = 1
    for degree_j in range(0, len(polyfit_roots)):
        # x[0]**n * p[0] + ... + x[0] * p[n-1] + p[n] = y[0]
        ygpz_val *= (jw - polyfit_roots[degree_j])
    # ygpz_amp_val = abs(ygpz_val)

    return yfit_val, ygpz_val


def fit_freq_response(plotphase, order, infile):
    order = int(order)
    arg = "[fit_freq_response] plotphase %s .." % plotphase
    print(arg)
    arg = "[fit_freq_response] order %s .." % order
    print(arg)

    # get FFT
    arg = "[fit_freq_response] ***************************************"
    print(arg)
    arg = "[fit_freq_response] getting FFT .."
    print(arg)
    nbits = 24
    dbfs_amplitude = (2.0**(nbits-1))
    freq, norm_abs_rfft_dbfs, norm_abs_rfft = ss_fft.do_rfft(infile=infile,
                                                             outfile=infile.replace(".MSEED", "_fft.png"),
                                                             nbits=nbits)
    freq = freq[:-2500]
    norm_abs_rfft_dbfs = norm_abs_rfft_dbfs[:-2500]
    norm_abs_rfft = norm_abs_rfft[:-2500]
    arg = "[fit_freq_response] evaluating only up to freq[len(freq)-1] %s" % freq[len(freq)-1]
    print(arg)

    # get polyfit
    arg = "[fit_freq_response] ***************************************"
    print(arg)
    arg = "[fit_freq_response] getting polyfit .."
    print(arg)
    polyfit_coeffs, residuals, rank, singular_values, rcond = np.polyfit(freq, norm_abs_rfft, order, full=True)
    polyfit_degrees = len(polyfit_coeffs) - 1
    arg = "[fit_freq_response] polyfit_coeffs %s" % polyfit_coeffs
    print(arg)
    sqrt_mse = np.sqrt(residuals[0]/len(freq))
    arg = "[fit_freq_response] sqrt(MSE) %s [Counts, not in dB scale]" % sqrt_mse
    print(arg)
    sqrt_mse_dbfs = 20*np.log10(sqrt_mse/(2.0**(nbits-1)))

    # get estimated transfer function in GainPolesZeros form
    polyfit_roots = np.roots(polyfit_coeffs)
    polyfit_roots = polyfit_roots*(2.0*math.pi)*1j    # convert them to rads/s form
    arg = "[fit_freq_response] polyfit_roots %s" % polyfit_roots
    print(arg)
    polyfit_poly = np.poly(polyfit_roots)
    arg = "[fit_freq_response] polyfit_poly %s" % polyfit_poly
    print(arg)
    arg = "[fit_freq_response] polyfit_coeff %s" % polyfit_coeffs
    print(arg)

    arg = "[fit_freq_response] ***************************************"
    print(arg)
    arg = "[fit_freq_response] Calculate Gain to fit input .."
    print(arg)
    # Calculate Gain to fit points
    # Evaluate Tranfer Function
    freq_tested = 0    # at 0 [Hz]
    yfit_val, ygpz_val = evaluate_functions(freq_tested, polyfit_roots, polyfit_degrees, polyfit_coeffs)
    yfit_amp_val = abs(yfit_val)
    ygpz_amp_val = abs(ygpz_val)
    # Calculate Gain
    ygpz_gain = yfit_amp_val/ygpz_amp_val
    arg = "[fit_freq_response] yfit_amp_val %s" % yfit_amp_val
    print(arg)
    arg = "[fit_freq_response] ygpz_amp_val %s" % ygpz_amp_val
    print(arg)
    arg = "[fit_freq_response] ygpz_gain %s" % ygpz_gain
    print(arg)
    # Show transfer function
    ygpz_transfer_function = ""
    for degree_j in range(0, len(polyfit_roots)):
        # x[0]**n * p[0] + ... + x[0] * p[n-1] + p[n] = y[0]
        ygpz_transfer_function += "(jw - %s)*" % polyfit_roots[degree_j]
    ygpz_transfer_function = "H(jw) = %s*%s" % (ygpz_gain, ygpz_transfer_function[:-1])
    arg = "[fit_freq_response] Scaled to fit input => %s" % ygpz_transfer_function
    print(arg)

    arg = "[fit_freq_response] ***************************************"
    print(arg)
    arg = "[fit_freq_response] Calculate Gain at 1 [Hz] .."
    print(arg)
    # Calculate Gain normalized 1 [Hz]
    # Evaluate Tranfer Function
    freq_tested = 1.0    # at 1 [Hz]
    yfit_val, ygpz_val = evaluate_functions(freq_tested, polyfit_roots, polyfit_degrees, polyfit_coeffs)
    yfit_amp_val = abs(yfit_val)
    ygpz_amp_val = abs(ygpz_val)
    # Calculate Gain
    ygpz_gain_1_at1hz = 1.0/ygpz_amp_val
    arg = "[fit_freq_response] yfit_amp_val %s" % yfit_amp_val
    print(arg)
    arg = "[fit_freq_response] ygpz_amp_val %s" % ygpz_amp_val
    print(arg)
    arg = "[fit_freq_response] ygpz_gain_1_at1hz %s" % ygpz_gain_1_at1hz
    print(arg)
    # Show transfer function
    ygpz_transfer_function = ""
    for degree_j in range(0, len(polyfit_roots)):
        # x[0]**n * p[0] + ... + x[0] * p[n-1] + p[n] = y[0]
        ygpz_transfer_function += "(jw - %s)*" % polyfit_roots[degree_j]
    ygpz_transfer_function = "H(jw) = %s*%s" % (ygpz_gain_1_at1hz, ygpz_transfer_function[:-1])
    arg = "[fit_freq_response] Normalized to 1.0 [Hz] => %s" % ygpz_transfer_function
    print(arg)

    # Evaluate estimated transfer function
    arg = "[fit_freq_response] ***************************************"
    print(arg)
    arg = "[fit_freq_response] Evaluating transfer function .."
    print(arg)
    arg = "[fit_freq_response] freq[0] = %s" % freq[0]
    print(arg)
    arg = "[fit_freq_response] freq[len(freq)-1] = %s" % freq[len(freq)-1]
    print(arg)
    yfit_amp_dbfs = []
    ygpz_amp_dbfs = []
    yfit_amp_db = []
    ygpz_amp_db = []
    yfit_phase = []
    ygpz_phase = []
    # freq_log10 = []
    for freq_i in freq:
        # Evaluate Transfer FUnctions at freq_i
        yfit_val, ygpz_val = evaluate_functions(freq_i, polyfit_roots, polyfit_degrees, polyfit_coeffs)
        yfit_amp_val = abs(yfit_val)
        ygpz_amp_val = abs(ygpz_val)
        yfit_phase_val = cmath.phase(yfit_val)
        ygpz_phase_val = cmath.phase(ygpz_val)

        # Apply calculated gain
        # ygpz_amp_val *= ygpz_gain                        # Count/Volt in dBFS scale, gain to fit to data
        ygpz_amp_val *= ygpz_gain_1_at1hz * (2**24/4.0)  # Count/Volt in dBFS scale, gain_1V_at1Hz

        # Create arrays to plot Amp-Freq and Phase-Freq graphs
        # Amplitude
        ygpz_amp_dbfs.append(20*np.log10(ygpz_amp_val/dbfs_amplitude))
        yfit_amp_dbfs.append(20*np.log10(yfit_amp_val/dbfs_amplitude))
        ygpz_amp_db.append(20*np.log10(ygpz_amp_val))
        yfit_amp_db.append(20*np.log10(yfit_amp_val))
        # yfit_amp_dbfs.append(yfit_amp_val)

        # Pahse
        y_phase_val = ygpz_phase_val*180.0/math.pi
        ygpz_phase.append(y_phase_val)
        y_phase_val = yfit_phase_val*180.0/math.pi
        yfit_phase.append(y_phase_val)

        # # Frequency axis
        # freq_log10 = 10*np.log10(freq)

    # arg = "[fit_freq_response] len(yfit_amp_dbfs) %s" % len(yfit_amp_dbfs)
    # print(arg)
    arg = "[fit_freq_response] yfit_amp_dbfs[0] %s" % yfit_amp_dbfs[0]
    print(arg)
    arg = "[fit_freq_response] yfit_amp_dbfs[len(yfit_amp_dbfs)-1] %s" % yfit_amp_dbfs[len(yfit_amp_dbfs)-1]
    print(arg)
    # arg = "[fit_freq_response] len(ygpz_amp_dbfs) %s" % len(ygpz_amp_dbfs)
    # print(arg)
    arg = "[fit_freq_response] ygpz_amp_dbfs[0] %s" % ygpz_amp_dbfs[0]
    print(arg)
    arg = "[fit_freq_response] ygpz_amp_dbfs[len(ygpz_amp_dbfs)-1] %s" % ygpz_amp_dbfs[len(ygpz_amp_dbfs)-1]
    print(arg)

    # Plot Amp-Freq and Phase-Freq graphs
    arg = "[fit_freq_response] ***************************************"
    print(arg)
    arg = "[fit_freq_response] getting plots .."
    print(arg)

    plt.clf()
    if plotphase:
        plt.subplot(2, 1, 1)

    # Plot Title
    plot_title = 'Polyfit of %sth order \n sqrt(MSE) = %s [Counts] => %s [dBFS]' \
                 % (polyfit_degrees, sqrt_mse, sqrt_mse_dbfs)

    # Plot amplitude
    # label1 = "System response"  # label1 = ygpz_transfer_function
    # plt.plot(freq, norm_abs_rfft_dbfs, marker='o', linestyle="", markersize=2, color='red', label=label1)
    # label_db = "Estimated Transfer Function dB"
    # plt.plot(freq, ygpz_amp_db, marker='o', color='yellow', markersize=3, label=label_db)
    ind = np.argmax(freq > 0.9999)
    val = freq[ind]
    val2 = ygpz_amp_dbfs[ind]
    arg = "[fit_freq_response] freq[%s] = %s => ygpz_amp_dbfs[%s] = %s" % (ind, val, ind, val2)
    print(arg)
    label_dbfs = "Estimated Transfer Function \n %s [Hz] => %s [dBFS]" % (val, val2)
    plt.plot(freq, ygpz_amp_dbfs, marker='o', color='blue', markersize=3, label=label_dbfs)
    # label3 = "Polyfit"
    # plt.plot(freq, yfit_amp_dbfs, marker='o', color='yellow', label=label3)
    plt.title(plot_title)
    plt.xlabel('Frequency [Hz]')
    # plt.xscale('log', basey=10, nonposy='clip')
    plt.ylabel('Amplitude Response: Counts / Volts [dBFS], %s bits' % nbits)
    # plt.ylabel('Amplitude Response: Counts / Volts [dB]')

    # Leyend
    # blue_line = mlines.Line2D([], [], color='blue', marker='*', markersize=15, label='Blue stars')
    # red_patch = mpatches.Patch(color='red', label='The red data')
    # plt.legend(handles=[red_patch])
    # plt.legend(handles=[blue_line])
    # plt.legend()
    # plt.legend(loc=2,prop={'size': 6})
    plt.legend(prop={'size': 10})

    # Extra info (max Freq, dBFS,  etc)
    # y_axis = np.float64(ygpz_amp_dbfs)
    # x_axis = np.float64(frequencies)
    # y_axis_argmax = int(freq_cutoff/freq_step) - 1
    # y_axis_max = y_axis[y_axis_argmax]
    # y_axis_min = y_axis.min()
    # x_axis_min = x_axis.min()
    # print('y_axis_argmax = %s' % y_axis_argmax)
    # print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
    # print('y_axis_max() = %s' % y_axis_max)
    # arg_freqinfo = '%s [Hz]\n %s [dB]' % (x_axis[y_axis_argmax], y_axis_max)
    # plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]-0.5, y_axis_max-0.5),
    #              xytext=(x_axis_min+1, y_axis_min+1),
    #              arrowprops=dict(facecolor='black', shrink=0.05),
    #              )

    # Plot phase
    if plotphase:
        plt.subplot(2, 1, 2)
        # ax[1].plot(frequencies_log10, ygpz_phase, marker='o')
        plt.plot(freq, ygpz_phase, marker='o', markersize=2)
        plt.xlabel('Frequency [Hz]')
        # plt.xscale('log', basey=10, nonposy='clip')
        plt.ylabel('Phase Response [degrees]')

        # Extra info (max Freq, dBFS,  etc)
        # y_axis = ygpz_phase
        # x_axis = frequencies
        # y_axis_argmax = int(freq_cutoff/freq_step) - 1
        # y_axis_max = y_axis[y_axis_argmax]
        # print('y_axis_argmax = %s' % y_axis_argmax)
        # print('y_axis[%s] = %s' % (y_axis_argmax, y_axis[y_axis_argmax]))
        # print('y_axis_max() = %s' % y_axis_max)
        # arg_freqinfo = ' %s [Hz]\n %s [degrees]' % (x_axis[y_axis_argmax], y_axis_max)
        # plt.annotate(arg_freqinfo, xy=(x_axis[y_axis_argmax]+1, y_axis_max+1),
        #              xytext=(x_axis[y_axis_argmax]+5, y_axis_max+10),
        #              arrowprops=dict(facecolor='black', shrink=0.05),
        #              )

    outfile = infile.replace(".MSEED", "_polyfit.png")
    arg = "[fit_freq_response] saving file %s .." % outfile
    print(arg)
    plt.savefig(outfile)
    # plt.show()
    plt.clf()

    return freq, ygpz_amp_dbfs, ygpz_phase, polyfit_roots, ygpz_gain

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot given file(s) (obspy wrapper)')
    parser.add_argument('--infile', action='store', help='name of input file')
    parser.add_argument('--order', action='store', help='function order')
    parser.add_argument('--plotphase', action='store_true', help='signal to generate')
    #parser.add_argument('file', type=argparse.FileType('r'), nargs='+')
    # parser.add_argument('--day_plot', action='store_true', help='day plot of the given file(s), normally same channel')
    args = parser.parse_args()

    if args.order is None:
        args.order = 4
    fit_freq_response(plotphase=args.plotphase, order=args.order, infile=args.infile)