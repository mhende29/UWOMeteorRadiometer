

from __future__ import print_function, division, absolute_import

import argparse
import os
import sys
import math
import csv
from datetime import datetime
from GetRDMConfig import RDMConfig, readConfig

import numpy as np
import scipy.signal
import scipy.optimize
import scipy.interpolate

from getRDMData import getRDMData

def initNotchFilter(sps, freq, band, ripple, order, filter_type):
    """ Initializes a noth filter with the given parameters.

    Arguments:
        sps: [float] Samples per second.
        freq: [float] Middle frequency of the filter in Hz.
        band: [float] Width of the band in Hz.
        ripple: [float] Ripple in the bandpass.
        order: [int] Filter order.
        filter_type: [str] Filter type, e.g. 'cheby1'

    Return:
        (b, a): [tuple] Filter parameters.

    """

    # Compute Nyquist frequency
    nyq  = sps/2.0

    # Lower frequency cutoff
    low  = freq - band/2.0

    # Upper frequency cutoff
    high = freq + band/2.0

    low  = low/nyq
    high = high/nyq

    # Init the filter
    b, a = scipy.signal.iirfilter(order, [low, high], rp=ripple, btype='bandstop', analog=False, 
        ftype=filter_type)

    return b, a



def butterworthBandpassFilter(lowcut, highcut, fs, order=5):
    """ Butterworth bandpass filter.

    Argument:
        lowcut: [float] Lower bandpass frequency (Hz).
        highcut: [float] Upper bandpass frequency (Hz).
        fs: [float] Sampling rate (Hz).

    Keyword arguments:
        order: [int] Butterworth filter order.

    Return:
        (b, a): [tuple] Butterworth filter.

    """

    # Calculate the Nyquist frequency
    nyq = 0.5*fs

    low = lowcut/nyq
    high = highcut/nyq

    # Init the filter
    b, a = scipy.signal.butter(order, [low, high], analog=False, btype='band')

    # Check if the filter is unstable
    if np.all(np.abs(np.roots(a)) < 1):
        print("""The filter with bands from {:.2f} Hz to {:.2f} Hz is unstable and should not be \
used! It's roots are smaller than 1.""".format(
            lowcut, highcut))

    return b, a



def filterBandpass(data, sps, bandpass_low, bandpass_high, order=6):
    """ Run bandpass filtering. """

    # Init the butterworth bandpass filter
    butter_b, butter_a = butterworthBandpassFilter(bandpass_low, bandpass_high, sps, order=order)

    # Filter the data
    waveform_data = scipy.signal.lfilter(butter_b, butter_a, np.copy(data))

    return waveform_data



def filterLP(data, sps, mains_freq, lowpass=True, filter_order=3, additional=None, low_pass_cutoff=500):
    """ Filter out the light pollution using Chebyshev filters. 
    
    Arguments:
        data: [ndarray] Unfiltered data.
        sps: [float] Samples per second.
        mains_freq: [float] Electric grid frequency in Hz (50 for Europe, 60 for NA).
    
    Keyword arguments:
        lowpass: [bool] Apply at lowpass filter (cutoff freq at the mains freq). True by default.
        filter_oder: [int] Order of the Chebyasev filter (3 by default).
        additional: [list] A list of (frequency, band width) tuples which define additional frequencies for
            filtering.

    Return:
        [ndarray]: Filtered data.

    """

    filter_type = 'cheby1'
    ripple = 10.0

    # Generate filter parameters for all harmonics
    filters_params = []
    for i in range(int((sps/2)/mains_freq)):

        # Compute the current harmonic frequency
        f_har = mains_freq*(i + 1)

        # If the lowpass filter is on, skip all frequencies above 5x the mains frequency
        if lowpass:
            if f_har > 2*low_pass_cutoff:
                continue

        # Set proper filter band width, depending on the harmonic number (first ones are wider)
        if i == 0:
            band_har = 7.5
        elif i == 1:
            band_har = 5.5
        else:
            band_har = 3.5

        filters_params.append([sps, f_har, band_har, ripple, filter_order, filter_type])


    if additional is not None:
        for freq, band in additional:

            # If the lowpass filter is on, skip all frequencies above 5x the mains frequency
            if lowpass:
                if freq > 2*low_pass_cutoff:
                    continue

            filters_params.append([sps, freq, band, ripple, filter_order, filter_type])

    filtered_data = np.copy(data)

    # Detrend the data
    filtered_data = scipy.signal.detrend(filtered_data)


    # Filter data using notch filters
    for filt_param in filters_params:

        print(filt_param)

        # Init the filter
        b, a = initNotchFilter(*filt_param)

        # Filter the data
        filtered_data = scipy.signal.lfilter(b, a, filtered_data)



    if lowpass:

        ### Apply a lowpass filter which will filter out everything above the grid frequency ###

        # Init the lowpass filter
        Wn = low_pass_cutoff/(sps/2.0)
        b, a = scipy.signal.butter(6, Wn)

        # Filter the data
        filtered_data = scipy.signal.filtfilt(b, a, filtered_data)

        ##############################


    return filtered_data


def movingAverage(arr, n=3):
    """ Perform a moving average on an array with the window size n.

    Arguments:
        arr: [ndarray] Numpy array of values.

    Keyword arguments:
        n: [int] Averaging window.

    Return:
        [ndarray] Averaged array. The size of the array is always by n-1 smaller than the input array.

    """

    ret = np.cumsum(arr, dtype=float)

    ret[n:] = ret[n:] - ret[:-n]

    return ret[n - 1:]/n




def datestr2UnixTime(time_str, UT_corr=0.0):
    """ Convert date and time to Unix time. 
    Arguments:
        time_str: [str]

    Kwargs:
        millisecond: [int] milliseconds (optional)
        UT_corr: [float] UT correction in hours (difference from local time to UT)
    
    Return:
        [float] Unix time

    """

    # Convert the time string to datetime
    dt = datetime.datetime.strptime(time_str, "%Y%m%d-%H%M%S.%f") - datetime.timedelta(hours=UT_corr)

    # UTC unix timestamp
    unix_timestamp = (dt - datetime.datetime(1970, 1, 1)).total_seconds()

    return unix_timestamp




def sine(t, a, f, phi, c):
    return np.abs(a)*np.sin(2*np.pi*f*t + phi%(2*np.pi)) + c



def fitSine(time_data, intensity_data, f0, p0 = None):
    """ Fits a sine to the given data.

    Arguments:
        f0: [float] Initial guess of the frequency.

    """ 

    a0 = np.std(intensity_data)
    phi0 = 0
    c0 = np.mean(intensity_data)

    # Initial guess
    if(p0 is None):
        p0 = [a0, f0, phi0, c0]

    popt, _ = scipy.optimize.curve_fit(sine, time_data, intensity_data, p0=p0)

    a, f, phi, c = popt
    popt[0] = np.abs(a)
    popt[2] = phi%(2*np.pi)

    return popt



def sineSlide(time_data, intensity_data, f0, window_width, shift_width):
    """ Fits a sines over windows that cover the entire signal to estimate the frequency drift.

    Arguments:
        time_data: [list of floats] The time stamp of the samples.
        intensity_data: [list of ints] The intensity of the data.
        f0: [float] Initial guess of the frequency.
        window_width: [float] Window width in seconds.
        shift_width: [float] How far over the window shifts for the next fitting.
    """ 

    times = []
    amplitudes = []
    frequencies = []
    phases = []
    offsets = []
    residuals = []

    time_relative = time_data - np.min(time_data)
    sps = len(time_relative)/(time_relative[-1] - time_relative[0])  
    width = math.ceil(sps*window_width)
    shift = math.ceil(sps*shift_width)

    loops = int((len(intensity_data) - width)//shift)

    print("Data points:", len(intensity_data))
    print("SPS:", sps)
    print("Samples per windows:", width)
    print("Samples per window shift:", shift)
    print("Number of window shifts:", loops)

    for i in range(loops):

        temp_time = time_relative[int(i*shift): int((i*shift) + width)]
        temp_data = intensity_data[int(i*shift): int((i*shift) + width)]

        if(i == 0):
            sine_fit = fitSine(temp_time, temp_data, f0=f0)
        else:
            sine_fit = fitSine(temp_time, temp_data, f0=f0, p0 = [amplitudes[i-1], frequencies[i-1], phases[i-1], offsets[i-1]])
        
        times.append([temp_time[0] + np.min(time_data), temp_time[-1] + np.min(time_data)])
        amplitudes.append(sine_fit[0])
        frequencies.append(sine_fit[1])
        phases.append(sine_fit[2])
        offsets.append(sine_fit[3])

        # Compute the sttdev of residuals
        temp_fit = list(sine_fit)
        temp_fit[-1] = 0
        res_data = temp_data - sine(temp_time, *sine_fit)

        # plt.plot(temp_time, temp_data, label='original')
        # plt.plot(temp_time, temp_data - sine(temp_time, *temp_fit), label='filtered')
        # plt.plot(temp_time, sine(temp_time, *sine_fit), label='fit')
        # plt.legend()
        # plt.show()

        residuals.append(np.std(res_data))

        if(i != (loops - 1)):
            print("Sines fitted: {:.2%}".format(i/(loops - 1)),end = "\r")
        else:
            print("Sines fitted: {:.2%}".format(i/(loops - 1)),end = "\n")

    # Unwrap the phase
    phases = np.unwrap(phases)

    coefs = [times, amplitudes, frequencies,phases, offsets, residuals]

    return coefs



def filterInterpolatedSines(time_data, intensity_data, sine_fits):


    sine_times_unix, amplitudes, frequencies, phases, offsets = sine_fits

    amp_interp = scipy.interpolate.PchipInterpolator(sine_times_unix, amplitudes)
    freq_interp = scipy.interpolate.PchipInterpolator(sine_times_unix, frequencies)
    phase_interp = scipy.interpolate.PchipInterpolator(sine_times_unix, phases)
    offset_interp = scipy.interpolate.PchipInterpolator(sine_times_unix, offsets)

    fitted_sines = sine(time_data, amp_interp(time_data), freq_interp(time_data), phase_interp(time_data), offset_interp(time_data))

    filtered_data = intensity_data - sine(time_data, amp_interp(time_data), freq_interp(time_data), phase_interp(time_data), np.zeros_like(time_data))


    return filtered_data, fitted_sines




def detrend(time_data, intensity_data, coefs):

    time_relative = time_data - np.min(time_data)

    for i in range(len(times)):

        sine_fit = [amplitudes[i], frequencies[i], phases[i], offsets[i]]

        lower_bound = times[i][0]
        upper_bound = times[i][0]

        # Find indices to cut
        beg_ind = (np.abs(np.array(time_data) - times[i][0])).argmin()
        end_ind = (np.abs(np.array(time_data) - times[i][1])).argmin()


        #sines = sine_fit['fitfunc'](time_relative[beg_ind:end_ind])

    #print(sines)




    return
    


def test_noise_removal(time_data, intensity_data):

    time_relative = time_data - np.min(time_data)
    sps = len(time_relative)/(time_relative[-1] - time_relative[0])
    time_chunk = 1/60
    samples_per_60Hz = int(math.ceil(sps*time_chunk))
    print(samples_per_60Hz)

    done = False
    noise_chunks = []
    loop_counter = 0

    while done is not True:
        current_chunk = intensity_data[loop_counter*samples_per_60Hz : (loop_counter + 1)*samples_per_60Hz]

        if(len(current_chunk) == 0):
            break

        current_chunk_array = np.array(current_chunk)

        if (((current_chunk_array.max() - current_chunk_array.min()) < 8000) and (len(current_chunk) == samples_per_60Hz)):
            noise_chunks.append(current_chunk)

        loop_counter += 1

    noise_values = np.matrix(noise_chunks)
    noise_avg = noise_values.mean(0) - (noise_values.mean(0)).min()
    noise_avg = [noise_avg.item(i) for i in range(36)]
    noise_avg = np.array(noise_avg)

    loop_counter = 0
    clean_data = []

    while True:

        if(((loop_counter + 1)*samples_per_60Hz) < (len(intensity_data)-1)):
            current_chunk = intensity_data[loop_counter*samples_per_60Hz : (loop_counter + 1)*samples_per_60Hz]
        else:
            current_chunk = intensity_data[loop_counter*samples_per_60Hz :]

        if(len(current_chunk) == 0):
            break

        cleaned = current_chunk - noise_avg[0:len(current_chunk)]
        clean_data.append(cleaned)


        loop_counter += 1

    #good_data = [data for data in *clean_data]
    clean_data  = np.concatenate(clean_data, axis=0)

    print(len(intensity_data))
    print(len(clean_data))

    return clean_data
 
def export_To_CSV(dir_path, station_code, station_channel, time_data, intensity_data, filtered):
    beg_time = datetime.utcfromtimestamp(time_data[0])
    end_time = datetime.utcfromtimestamp(time_data[-1])

    waiting_time = end_time - beg_time
    time_dif = int(waiting_time.total_seconds())
    
    if(filtered):
        file_name = "{:s}_{:s}_{:04d}{:02d}{:02d}-{:02d}{:02d}{:02d}.{:06d}_{:02d}{:02d}{:02d}.{:06d}_filtered.csv".format(station_code, station_channel, beg_time.year,beg_time.month,beg_time.day,beg_time.hour,beg_time.minute,beg_time.second, beg_time.microsecond, end_time.hour,end_time.minute,end_time.second, end_time.microsecond)
    else:
        file_name = "{:s}_{:s}_{:04d}{:02d}{:02d}-{:02d}{:02d}{:02d}.{:06d}_{:02d}{:02d}{:02d}.{:06d}.csv".format(station_code, station_channel, beg_time.year,beg_time.month,beg_time.day,beg_time.hour,beg_time.minute,beg_time.second, beg_time.microsecond, end_time.hour,end_time.minute,end_time.second, end_time.microsecond)
    
    header = ["# Unix Times", "Intensity"]

    data = [[time_data[i], intensity_data[i]] for i in range(len(time_data))]

    with open(os.path.join(dir_path, file_name), 'w') as csvfile:
        file = csv.writer(csvfile)

        file.writerow(header)

        file.writerows(data)


if __name__ == "__main__":
        
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    
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

    arg_p.add_argument('-e', action = 'store_true', help="""If enabled produces and exports a csv file.""")

    arg_p.add_argument('-n', action = 'store_true', help="""Loads the night plots.""")
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

    if(cml_args.e is True):
        if(not os.path.isdir(csv_path)):
            os.mkdir(csv_path, 0o755)
        export_To_CSV(csv_path, cml_args.code, cml_args.channel, unix_times, intensity, filtered)
        filtered = True

    if(cml_args.n is True):
        yr_mon_day = cml_args.time[:8]

        night = [folder for folder in os.listdir(archived_data_path) if(folder.startswith(cml_args.code  + "_" + cml_args.channel + "_" + yr_mon_day))]

        pngs = [plot for plot in os.listdir(os.path.join(archived_data_path, *night)) if(plot.endswith(".png"))]

        night_plot = os.path.join(os.path.join(archived_data_path, *night),pngs[0])
        max_minus = os.path.join(os.path.join(archived_data_path, *night),pngs[1])

        img = mpimg.imread(night_plot)
        fig = plt.imshow(img)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        plt.tight_layout()
        plt.box(on=None)
        plt.show()

        img = mpimg.imread(max_minus)
        fig = plt.imshow(img)
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        plt.tight_layout()
        plt.box(on=None)
        plt.show()

    # Compute relative time since the beginning of the recording
    time_relative = unix_times - np.min(unix_times)
    
    # Create a list containing all times as datetime objects
    all_datetime = [datetime.utcfromtimestamp(t) for t in unix_times]

    # Compute samples per second
    sps = len(time_relative)/(time_relative[-1] - time_relative[0])
    print('SPS:', sps)

    #clean_data = test_noise_removal(unix_times, intensity)

    #plt.plot(unix_times, intensity)
    #plt.plot(unix_times, clean_data)
    #plt.show()

    
    # filtered_data = np.array(intensity)

    # for i in range(1):

    #     # Fit the sines by sliding a window
    #     coefs = sineSlide(unix_times, filtered_data, f0 = (i + 1)*60, window_width = 1.0, shift_width = 0.1)
        

    #     #detrend(unix_times, intensity, coefs)


    #     times, amplitudes, frequencies,phases, offsets, residuals = coefs

    #     sine_times_unix = [t for t in np.mean(np.array(times), axis = 1)]
    #     sine_times = [datetime.utcfromtimestamp(t) for t in sine_times_unix]

    #     sine_fits = [sine_times_unix, amplitudes, frequencies, phases, offsets]

    #     # Filter original data with fitted sines
    #     filtered_data, fitted_sines = filterInterpolatedSines(unix_times, filtered_data, sine_fits)


    # plt.plot(unix_times, intensity, linewidth=0.2)
    # plt.plot(unix_times, filtered_data, linewidth=0.2)
    # #plt.plot(unix_times, fitted_sines,linewidth=0.2)
    # plt.show()


    # fig, (ax1, ax2) = plt.subplots(nrows=2)
    # ax1.specgram(intensity, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=16)
    # ax2.specgram(filtered_data, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=16)

    # plt.show()



    # #######################################################################################################################################################################################



    # plt.plot(sine_times, frequencies, linewidth=0.2)
    
    # if(cml_args.time.find(".") is not -1):
    #     title = ("Frequency profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S.%f")).strftime("%H:%M:%S.%f on %B%e, %Y UTC")))
    #     micros = title[title.find(".") + 1:title.find(".") + 6 + 1]
    #     newmicros = micros.rstrip("0")
    #     title = title.replace(micros, newmicros)
    # else:
    #     title = "Frequency profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S")).strftime("%H:%M:%S on %B%e, %Y UTC"))

    # plt.title(title)
    # plt.ylabel('Frequency')
    # plt.xlabel('Time')
    # plt.xticks(rotation=30)
    # plt.show()

    # plt.plot(sine_times, amplitudes, linewidth=0.2)
    
    # if(cml_args.time.find(".") is not -1):
    #     title = ("Amplitude profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S.%f")).strftime("%H:%M:%S.%f on %B%e, %Y UTC")))
    #     micros = title[title.find(".") + 1:title.find(".") + 6 + 1]
    #     newmicros = micros.rstrip("0")
    #     title = title.replace(micros, newmicros)
    # else:
    #     title = "Amplitude profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S")).strftime("%H:%M:%S on %B%e, %Y UTC"))

    # plt.title(title)
    # plt.ylabel('Amplitude')
    # plt.xlabel('Time')
    # plt.xticks(rotation=30)
    # plt.show()

    # plt.plot(sine_times, phases, linewidth=0.2)
    
    # if(cml_args.time.find(".") is not -1):
    #     title = ("Phase profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S.%f")).strftime("%H:%M:%S.%f on %B%e, %Y UTC")))
    #     micros = title[title.find(".") + 1:title.find(".") + 6 + 1]
    #     newmicros = micros.rstrip("0")
    #     title = title.replace(micros, newmicros)
    # else:
    #     title = "Phase profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S")).strftime("%H:%M:%S on %B%e, %Y UTC"))

    # plt.title(title)
    # plt.ylabel('Phase')
    # plt.xlabel('Time')
    # plt.xticks(rotation=30)
    # plt.show()

    # plt.plot(sine_times, offsets, linewidth=0.2)
    
    # if(cml_args.time.find(".") is not -1):
    #     title = ("Offset profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S.%f")).strftime("%H:%M:%S.%f on %B%e, %Y UTC")))
    #     micros = title[title.find(".") + 1:title.find(".") + 6 + 1]
    #     newmicros = micros.rstrip("0")
    #     title = title.replace(micros, newmicros)
    # else:
    #     title = "Offset profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S")).strftime("%H:%M:%S on %B%e, %Y UTC"))

    # plt.title(title)
    # plt.ylabel('Offset')
    # plt.xlabel('Time')
    # plt.xticks(rotation=30)
    # plt.show()



    # plt.plot(sine_times, residuals, linewidth=0.2)

    # if(cml_args.time.find(".") is not -1):
    #     title = ("Residuals profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S.%f")).strftime("%H:%M:%S.%f on %B%e, %Y UTC")))
    #     micros = title[title.find(".") + 1:title.find(".") + 6 + 1]
    #     newmicros = micros.rstrip("0")
    #     title = title.replace(micros, newmicros)
    # else:
    #     title = "Residuals profile within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S")).strftime("%H:%M:%S on %B%e, %Y UTC"))

    # plt.title(title)
    # plt.ylabel('Stddev')
    # plt.xlabel('Time')
    # plt.xticks(rotation=30)
    # plt.show()
    
    # sys.exit()
    #######################################################################################################################################################################################

    # Plot raw data
    # Checks if there's a period in the given string in order to decide if micro-seconds should be incorporated.
    if(cml_args.time.find(".") is not -1):
        title = "Unfiltered data within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S.%f")).strftime("%H:%M:%S.%f on %B%e, %Y UTC")) 
        micros = title[title.find(".") + 1:title.find(".") + 6 + 1]
        newmicros = micros.rstrip("0")
        title = title.replace(micros, newmicros)
    else:
        title = "Unfiltered data within {:0.6g}s of {:s}".format(float(cml_args.range)/2,(datetime.strptime(cml_args.time,"%Y%m%d-%H%M%S")).strftime("%H:%M:%S on %B%e, %Y UTC"))

    # Plot the spectrogram
    ax1 = plt.subplot(211)
    plt.specgram(intensity, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=0)
    plt.title(title)
    plt.ylabel('Frequency (Hz)')
    plt.tick_params(bottom=False, labelbottom=False)
    plt.subplots_adjust(hspace=0)

    ax2 = plt.subplot(212, sharex = ax1)
    plt.plot(time_relative, intensity,linewidth=0.2)
    plt.ylabel('ADU')
    plt.xlabel('Time')
    plt.xlim([time_relative[0],time_relative[-1]])
    plt.tick_params(labelbottom=False)
    ax3 = ax2.twiny()
    ax3.set_xlim(all_datetime[0],all_datetime[-1])
    ax3.xaxis.set_ticks_position('bottom')
    
    plt.xticks(rotation=30)
    
    
    # plt.savefig('/home/pi/Desktop/{:s}.png'.format(file_name), dpi=300)
    
    plt.show()

    # Plot power spectral density
    #plt.psd(intensity, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    #plt.show()


    # Filter the light pollution
    mains_freq = 60.0 # Hz
    filtered_data = filterLP(intensity, sps, mains_freq, additional=[(160, 2.0), (32, 2.0), (94, 2.0), (553.0, 2.0), (614, 2.0), (40.0, 2.0), (20.0, 2.0)], lowpass=True)

    # Skip the first 5% samples
    beg_cut = int(len(time_relative)*0.05)
    time_relative = time_relative[beg_cut:] - time_relative[beg_cut]
    filtered_data = filtered_data[beg_cut:]
    unix_times = unix_times[beg_cut:]

    # Skip the last 1% samples
    end_cut = len(time_relative) - int(len(time_relative)*0.01)
    time_relative = time_relative[:end_cut]
    filtered_data = filtered_data[:end_cut]
    unix_times = unix_times[:end_cut] 

    # Gather datetime equivalents of the unix times
    unix_datetimes = [datetime.utcfromtimestamp(t) for t in unix_times]



    # # Apply a moving average
    # window_size = 3
    # filtered_data = movingAverage(rdm.intensity.astype(np.float64), n=window_size)
    # time_relative = time_relative[:len(filtered_data)]
    # filtered_data = filtered_data[::window_size]
    # time_relative = time_relative[::window_size]
    # sps = sps/window_size


    # # Apply a broad lowpass and a highpass filter
    # bandpass_low = 5.0 # Hz
    # bandpass_high = 60.0
    # filtered_data = filterBandpass(rdm.intensity.astype(np.float64), sps, bandpass_low, bandpass_high, order=3)

    # Print the filtered data
    print(filtered_data)


    # Plot filtered spectrogram
    ax1 = plt.subplot(211)
    plt.specgram(filtered_data, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=0)
    title = title.replace("Unfiltered", "Filtered")
    plt.title(title)
    plt.tick_params(bottom=False, labelbottom=False)
    plt.ylabel('Frequency (Hz)')
    plt.subplots_adjust(hspace=0)
    # Plot power spectral density
    #plt.psd(filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    #plt.show()

    mean_std = (np.mean(filtered_data) + np.std(filtered_data))*np.ones_like(time_relative)

    # Plot filtered data
    ax2 = plt.subplot(212, sharex = ax1)
    plt.plot(time_relative, filtered_data, label='Filtered Data')
    plt.plot(time_relative, mean_std, linestyle='--', label=r'Two $\sigma$ threshold')
    plt.legend(loc = 0)
    plt.xlabel('Time (s)')
    plt.ylabel('ADU')
    plt.xlim([time_relative[0],time_relative[-1]])
    plt.tick_params(labelbottom=False)
    ax3 = ax2.twiny()
    ax3.set_xlim(unix_datetimes[0],unix_datetimes[-1])
    ax3.xaxis.set_ticks_position('bottom')
    plt.xticks(rotation=30)
    plt.show()

    if(cml_args.e is True):
        if(not os.path.isdir(csv_path)):
            os.mkdir(csv_path, 0o755)
        export_To_CSV(csv_path, cml_args.code, cml_args.channel, unix_times, intensity, filtered)