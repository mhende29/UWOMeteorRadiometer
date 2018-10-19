
import argparse
import os
import sys
import math
import csv

from datetime import datetime
from GetRDMConfig import RDMConfig, readConfig
from AnalyzeData import psd_fun, back_fun, export_To_CSV
from getRDMData import getRDMData

import numpy as np
import scipy.signal
import scipy.optimize

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

if __name__ == "__main__":

    archived_data_path = os.path.expanduser("~/RadiometerData/ArchivedData")
    work_dir = os.getcwd()
    home_dir = os.path.dirname(work_dir)
    csv_storage_dir = "Exported_Data"
    csv_path = os.path.join(home_dir,csv_storage_dir)

    # Set up input arguments
    arg_p = argparse.ArgumentParser(description="Get radiometer files.")

    arg_p.add_argument('code', metavar='CODE', nargs='?', \
        help='Time of the event in the YYYYMMDD-HHMMSS.ms format.', type=str, default=None)

    arg_p.add_argument('channel', metavar='CHANNEL', nargs='?', \
        help='Time of the event in the YYYYMMDD-HHMMSS.ms format.', type=str, default=None)

    arg_p.add_argument('time', metavar='TIME', nargs='?', \
        help='Time of the event in the YYYYMMDD-HHMMSS.ms format.', type=str, default=None)

    arg_p.add_argument('range', metavar='DURATION_SECONDS', help="""Grabs data so the total, 
        covers the range/2 on both sides of time. """)

    # Parse input arguments
    cml_args = arg_p.parse_args()

    # Check if there is a config file in the library dir
    if(os.path.isfile(os.path.join(work_dir, "config.txt"))):

        # Read the config in the lib path
        analysis_config = readConfig(os.path.join(work_dir, "config.txt"))

        # Check if the server flag is set in the config
        if(analysis_config.read_from_server):
            archived_data_path = os.path.join(os.path.join("/home", "rdm_" + cml_args.code.lower()), "files")
    
    
    if not os.path.exists(archived_data_path):
        print('The archived data path: {:s} does not exist!'.format(archived_data_path))

    # Gather the radiometric data and the time stamps around the given time period
    intensity, unix_times = getRDMData(archived_data_path, cml_args.code, cml_args.channel, cml_args.time, cml_args.range)

    filtered = False

    if(not os.path.isdir(csv_path)):
        os.mkdir(csv_path, 0o755)
    export_To_CSV(csv_path, cml_args.code, cml_args.channel, unix_times, intensity, filtered)
    filtered = True

    sys.exit()
    # Compute relative time since the beginning of the recording
    time_relative = unix_times - np.min(unix_times)
    
    # Create a list containing all times as datetime objects
    all_datetime = [datetime.utcfromtimestamp(t) for t in unix_times]

    # Compute samples per second
    sps = len(time_relative)/(time_relative[-1] - time_relative[0])
    print('SPS:', sps)

    # Design notch filter
    # Filtering setup

    # Redefine the sampling frequency for simplicity
    fs = sps  # Sample frequency (Hz)

    # Calculate the original power spectral density function and clear the figure
    p_xx, freqs = plt.psd(intensity, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    plt.clf()

    # Generate an array that contains the frequencies found in the mains hum
    # Mains hum being the 60 Hz interference and its higher order harmonics 
    mains_hum = np.arange(60,fs/2,60)

    # Convert the powers into decibels (dB)
    power_db = 10*np.log10(p_xx)

    # Create a list of indices where all those frequencies are
    index = []
    for i in range(len(mains_hum)):
        index = np.append(index, np.abs(freqs - mains_hum[i]).argmin()).astype(int)

    mean_indices = np.copy(index)
    mean_indices = np.append(np.array([0]), mean_indices)

    # Define a list to store the indices of the maximum powers, which are near but not exactly the same as the main frequency indices
    max_powers = []

    # Define a list to store the indices of the powers that aren't the mains, in order to make a fit that doesn't depend on the peaks
    fit_powers = list(np.copy(power_db))
    fit_freqs = list(np.copy(freqs))

    # Search 10 elements on either side of the nearest mains hum for the peak power index
    search_width = 10 
    for i in range(len(mains_hum)):
        max_powers = np.append(max_powers,index[i] - search_width + (p_xx[index[i]-search_width:index[i]+search_width]).argmax()).astype(int)

    for i in range(len(mains_hum)):
        j = len(mains_hum) - i - 1
        del fit_powers[max_powers[j] - search_width:max_powers[j] + search_width]
        del fit_freqs[max_powers[j] - search_width:max_powers[j] + search_width]

    # Begin a least-squares optimization fit that follows the noise in the psd
    # Begin initial guess's, 0 where additions/subtractions are and 1 where multiplications/divisions are
    x0 = np.array([0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0])
    
    # The optimization function in question is a three term gaussian model
    # f(x) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + d
    res_lsq = scipy.optimize.least_squares(psd_fun, x0, loss='cauchy', args=(fit_freqs, fit_powers))

    # Calculate the mean power using the least-squares solution
    power_means = psd_fun(res_lsq.x, freqs)

    # Calculate the standard deviation
    power_std = np.std(power_db)
    power_lower_lim = power_means - power_std/4.0
    power_upper_lim = power_means + power_std/4.0

    # Begin recursive filtering 

    filtered_data = np.copy(intensity)
    good_filtered_data = np.copy(filtered_data)

    i = 0
    for frequency in max_powers:
        power_to_verify = power_db[frequency]
        lower_bound = power_lower_lim[frequency]
        upper_bound = power_upper_lim[frequency]
        w0 = mains_hum[i]/(fs/2)
        Q = 1

        print("Filtering {:}Hz".format(int(mains_hum[i])))
        while((power_to_verify > upper_bound) or (power_to_verify < lower_bound)):
            b, a = scipy.signal.iirnotch(w0, Q)
            filtered_data = scipy.signal.lfilter(b, a, good_filtered_data)

            p_xx, f = plt.psd(filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
            #plt.show()
            plt.clf()
            power_db = 10*np.log10(p_xx)
            power_to_verify = power_db[frequency]

            if((power_to_verify < upper_bound) and (power_to_verify > lower_bound)):
                good_filtered_data = np.copy(filtered_data)
            else:
                filtered_data = np.copy(good_filtered_data)
                Q += mains_hum[i]/60
        i += 1

    corrected_unix_times = unix_times[int(0.02*len(unix_times)):-1]
    corrected_good_data = good_filtered_data[int(0.02*len(good_filtered_data)):-1]



    print("All frequencies filtered!")

    export_To_CSV(csv_path, cml_args.code, cml_args.channel, corrected_unix_times, corrected_good_data, filtered)

    sys.exit()