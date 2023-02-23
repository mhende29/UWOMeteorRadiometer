

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

import matplotlib.mlab

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


def generateEventName(station_code, station_channel, unix_times):
    """ Generate a template file name for the event. """

    beg_time = datetime.utcfromtimestamp(unix_times[0])
    end_time = datetime.utcfromtimestamp(unix_times[-1])

    event_name = "{:s}_{:s}_{:04d}{:02d}{:02d}-{:02d}{:02d}{:02d}.{:06d}_{:02d}{:02d}{:02d}.{:06d}".format(
        station_code, 
        station_channel, 
        beg_time.year,
        beg_time.month,
        beg_time.day,
        beg_time.hour,
        beg_time.minute,
        beg_time.second, 
        beg_time.microsecond, 
        end_time.hour,
        end_time.minute,
        end_time.second, 
        end_time.microsecond
        )

    return event_name

    
 
def exportCSV(dir_path, station_code, station_channel, unix_times, intensity_data, filtered):

    # Create an event name
    event_name = generateEventName(station_code, station_channel, unix_times)
    
    if(filtered):
        file_name = "{:s}_filtered.csv".format(event_name)
    else:
        file_name = "{:s}.csv".format(event_name)
    
    header = ["# Unix Times", "Intensity"]

    data = [[unix_times[i], intensity_data[i]] for i in range(len(unix_times))]

    with open(os.path.join(dir_path, file_name), 'w') as csvfile:
        file = csv.writer(csvfile)

        file.writerow(header)

        file.writerows(data)



def backgroundPSDModel(coefs, x, y=None):
    """  Three term gaussian model used to model the background frequency.
        f(x) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + d
        
    """

    # Compute the function value
    val = coefs[0] \
        + coefs[1]*np.exp(-((x - coefs[2])/coefs[3])**2) \
        + coefs[4]*np.exp(-((x - coefs[5])/coefs[6])**2) \
        + coefs[7]*np.exp(-((x - coefs[8])/coefs[9])**2)
    
    # If the data are given, compute the residuals
    if y is not None: 
        return val - y
    
    else:
        return val
    
def backgroundPSDModelResidualSum(coefs, x, y):

    # Squared value of each residual
    res = backgroundPSDModel(coefs, x, y=y)**2

    # Smooth approximation of l1 (absolute value) loss
    return np.sum(2*((1 + res)**0.5 - 1))


def linearBackgroundModel(coefs, x, y=None):

    # Compute the function value
    val = coefs[0] + coefs[1]*x

    # If the data are given, compute the residuals
    if y is not None: 
        return val - y
    
    else:
        return val



def showNightPlot(event_time, station_code, station_channel):
    """ Given the time and the station, show the night plot. """

    yr_mon_day = event_time[:8]

    # Find the appropriate night folder
    night = [folder for folder in os.listdir(archived_data_path) if folder.startswith(cml_args.code \
        + "_" + cml_args.channel + "_" + yr_mon_day)]

    if not night:
        print('The night was not found, night plot will not be shown!')
        return False

    # Get a list of PNGs in a given folder
    pngs = [plot for plot in os.listdir(os.path.join(archived_data_path, *night)) if(plot.endswith(".png"))]

    if not pngs:
        print('The night plot PNG was not found!')
        return False

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



if __name__ == "__main__":
        
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    
    archived_data_path = os.path.expanduser("~/RadiometerData/ArchivedData")

    work_dir = os.getcwd()
    home_dir = os.path.dirname(work_dir)
    csv_storage_dir = "Exported_Data"
    csv_path = os.path.join(home_dir,csv_storage_dir)

    # Create a directory on the disk for exporting events
    if(not os.path.isdir(csv_path)):
        os.mkdir(csv_path, 0o755)

    # Set up input arguments
    arg_p = argparse.ArgumentParser(description="Get radiometer files.")

    arg_p.add_argument('code', metavar='CODE', nargs='?', \
        help="""Radiometer station code, e.g. US0001. The first two letters are the ISO code of the country, 
        and the last 4 characters is the alphanumeric code of the station.""", type=str, default=None)

    arg_p.add_argument('channel', metavar='CHANNEL', nargs='?', \
        help='Channel of the station. Only single letters are accepted, such as A, B, C, etc.', type=str, \
        default=None)

    arg_p.add_argument('time', metavar='TIME', nargs='?', \
        help="""Time of the event in the YYYYMMDD-HHMMSS.ms format. Data will be centred around this time. 
        """, type=str, default=None)

    arg_p.add_argument('range', metavar='DURATION_SECONDS', type=float, nargs=1,
        help="""Duration in seconds of the data chunk that will be taken. The data will be taken in the range 
        of (-range/2, range/2), centered around the given time.""")

    arg_p.add_argument('-f', '--mainsfreq', metavar='MAINS_FREQ', nargs=1, \
        help="Frequency of the mains hum.", type=float)

    arg_p.add_argument('-e', '--exportcsv', action = 'store_true', \
        help="""If enabled produces and exports a csv file.""")

    arg_p.add_argument('-n', '--nightplot', action = 'store_true', \
        help="""Loads the night plots. They show an overview of all recorded intensities during the night.""")

    arg_p.add_argument('-p', '--showplots', action = 'store_true', \
        help="""Show the plots on the screen.""")
    
    arg_p.add_argument('-a', '--archivepath', metavar='ARCHIVE_PATH', type=str,
        help="""Path to the directory with archived. If not given, the default path will be used.""",
        default=None)
    
    arg_p.add_argument('-i', '--integrate', action = 'store_true', \
                       help="""If enabled, the program will attempt to automatically find the fireball and 
                       compute the integrated area under the curve.""")
    
    arg_p.add_argument('--addharms', type=float, nargs='+', \
        help="""List of additional harmonics to remove from the spectrum. Only the first frequency should be given, the higher order harmonics will be assumed to be multiples of the mains frequency.""")
    

    # Parse input arguments
    cml_args = arg_p.parse_args()


    cml_args.range = cml_args.range[0]


    ##########################################################################################################

    # Check if there is a config file in the library dir
    if(os.path.isfile(os.path.join(work_dir, "config.txt"))):

        # Read the config in the lib path
        config = readConfig(os.path.join(work_dir, "config.txt"))

        # Check if the server flag is set in the config
        if(config.read_from_server):
            archived_data_path = os.path.join(os.path.join("/home", "rdm_" + cml_args.code.lower()), "files")


    # If the config does not exist, load defualt values
    else:
        config = RDMConfig()



    # Assign the mains frequency hum if given
    if cml_args.mainsfreq:
        config.mains_frequency = cml_args.mainsfreq[0]


    # Take the archive path from the command line if given
    if cml_args.archivepath:
        archived_data_path = cml_args.archivepath
    
    
    if not os.path.exists(archived_data_path):
        print('The archived data path: {:s} does not exist!'.format(archived_data_path))


    # Extract the additional harmonics to remove from the spectrum
    additional_harmonics = None
    if cml_args.addharms:
        additional_harmonics = cml_args.addharms

    print("Additional harmonics to remove: ", additional_harmonics)

    # Gather the radiometric data and the time stamps around the given time period
    intensity, unix_times = getRDMData(archived_data_path, cml_args.code, cml_args.channel, cml_args.time, cml_args.range)

    filtered = False

    # Make a name for the event
    event_name = generateEventName(cml_args.code, cml_args.channel, unix_times)


    # Export data to CSV
    if cml_args.exportcsv:

        exportCSV(csv_path, cml_args.code, cml_args.channel, unix_times, intensity, filtered)

        filtered = True


    # Show the night plot
    if cml_args.nightplot:

        showNightPlot(cml_args.time, cml_args.code, cml_args.channel)



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
    # Mains hum being the 60/50 Hz interference and its higher order harmonics 
    mains_hum = np.arange(config.mains_frequency, fs/2, config.mains_frequency)

    # Add additional harmonics if given
    if additional_harmonics is not None:
        for harm in additional_harmonics:

            # Generate frequencies of the harmonics, with the same multiple as the mains hum
            harm_hum = np.arange(harm, fs/2, config.mains_frequency)

            # Add the harmonics to the list of mains hum frequencies
            mains_hum = np.append(mains_hum, harm_hum)

    # Sort the list of mains hum frequencies
    mains_hum = np.sort(mains_hum)


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
    # This will be identified as the center of the mains hum harmonic
    search_width = 10 
    for i in range(len(mains_hum)):
        max_powers = np.append(max_powers,index[i] - search_width + (p_xx[index[i]-search_width:index[i]+search_width]).argmax()).astype(int)

    for i in range(len(mains_hum)):
        j = len(mains_hum) - i - 1
        del fit_powers[max_powers[j] - search_width:max_powers[j] + search_width]
        del fit_freqs[max_powers[j] - search_width:max_powers[j] + search_width]

    

    # Begin a least-squares optimization fit that follows the noise in the psd
    x0 = np.array([np.median(fit_powers), # Background
                   np.mean(fit_powers), 2.0, 5.0, # First gaussian (mostly DC component)
                   0.0,   config.mains_frequency, 1.0, # Second Gaussian
                   0.0, 2*config.mains_frequency, 1.0  # Third Gaussian
                   ])
    

    # Fit the function using Nelder-Mead
    res = scipy.optimize.minimize(backgroundPSDModelResidualSum, x0, args=(fit_freqs, fit_powers), 
                                      method='Nelder-Mead')
    
    ### OLD CODE ###
    # Begin initial guess's, 0 where additions/subtractions are and 1 where multiplications/divisions are
    #x0 = np.array([0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0])
    
    # The optimization function in question is a three term gaussian model
    # f(x) = a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + a3*exp(-((x-b3)/c3)^2) + d

    # # Fit the function using least-squares
    # res = scipy.optimize.least_squares(backgroundPSDModel, x0, loss='soft_l1', args=(fit_freqs, fit_powers))

    ### ###
    
    print("Three term gaussian fit on the PSD of the background noise")
    print(res.x)

    # Calculate the mean power using the least-squares solution
    power_means = backgroundPSDModel(res.x, freqs)

    # Calculate the standard deviation of the fit
    power_std = np.std(power_db)
    power_lower_lim = power_means - power_std/4.0
    power_upper_lim = power_means + power_std/4.0

    # Plot the raw signal
    plt.plot(freqs, power_db, label='Raw signal PSD')

    # Plot the fit with a half-sigma sigma confidence bound for filtering 
    plt.plot(freqs, power_means, 'g', label='Background power fit')
    plt.plot(freqs, power_lower_lim, linestyle = '--', color = 'm', label=r'Confidence region: $\frac{1}{4}$$\sigma$')
    plt.plot(freqs, power_upper_lim, linestyle = '--', color = 'm')
    
    # Plot the identified peaks of the mains hum
    plt.plot(freqs[max_powers], 10*np.log10(p_xx[max_powers]), 'r+', label='Mains hum peaks')

    plt.title("Initial Power spectral Density of Raw Signal")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power Spectral Density (dB/Hz)")
    plt.grid(which='both')

    plt.legend(loc=0)

    plt.savefig(os.path.join(csv_path, event_name + "_psd.png"), dpi=150)

    if cml_args.showplots:
        plt.show()
    else:
        plt.clf()
        plt.close()

    
    
    # Begin recursive filtering of mains hum harmonics
    filtered_data = np.copy(intensity)
    good_filtered_data = np.copy(filtered_data)

    i = 0
    for peak_freq in max_powers:

        # Read the power at the peak
        power_to_verify = power_db[peak_freq]

        # Read the boundaries of the confidence region
        lower_bound = power_lower_lim[peak_freq]
        background_power = power_means[peak_freq]
        upper_bound = power_upper_lim[peak_freq]

        # Define the notch filter parameters
        w0 = mains_hum[i]/(fs/2)
        Q = 1.0

        print("Filtering {:}Hz".format(int(mains_hum[i])))


        def applyNotchFilter(data, w0, Q):
            """ Apply a notch filter to the data at the considered frequency
            
            Arguments:
                data: [ndarray] The data to be filtered
                w0: [float] The normalized frequency to be filtered
                Q: [float] The quality factor of the filter
                
            Returns:
                [ndarray] The filtered data

            """

            # Apply a notch filter to the data at the considered frequency
            b, a = scipy.signal.iirnotch(w0, Q)
            return scipy.signal.filtfilt(b, a, data)

        # Treat the quality factor as a variable to be optimized. Minimize the difference between the 
        # predicted background power and the filtered power at the peak frequency
        def filterResiduals(Q, w0, data, background_power, peak_freq):

            # Apply a notch filter to the data at the considered frequency
            filtered_data = applyNotchFilter(data, w0, Q)
            
            # Compute the PSD of the filtered data)
            #p_xx, f = plt.psd(filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
            #plt.clf()
            p_xx, _ = matplotlib.mlab.psd(filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)

            # Compute the power in dB
            power_db = 10*np.log10(p_xx)

            # Return the difference between the background power and the filtered power at the peak frequency
            return (background_power - power_db[peak_freq])**2

        # Begin the optimization
        res = scipy.optimize.minimize(filterResiduals, Q, 
                                      args=(w0, good_filtered_data, background_power, peak_freq),
                                      method='Nelder-Mead')
        
        print(" - Optimized notch filter Q = {:.2f}".format(res.x[0]))
        
        # Apply the optimized filter
        good_filtered_data = applyNotchFilter(good_filtered_data, w0, res.x[0])


        # ### OLD METHOD ###
        # # Increase the notch filter quality factor until the power at the peak is within the confidence region
        # while (power_to_verify > upper_bound) or (power_to_verify < lower_bound):

        #     # Apply a notch filter to the data at the considered frequency
        #     b, a = scipy.signal.iirnotch(w0, Q)
        #     #filtered_data = scipy.signal.lfilter(b, a, good_filtered_data) # old code, broken now

        #     filtered_data = scipy.signal.filtfilt(b, a, good_filtered_data)
            
        #     # Compute the PSD of the filtered data
        #     p_xx, f = plt.psd(filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
            
        #     # Read the power at the peak frequency
        #     power_db = 10*np.log10(p_xx)
        #     power_to_verify = power_db[peak_freq]

        #     # # Plot the power to verify and the confidence region
        #     # plt.plot(peak_freq, power_to_verify, 'r+')
        #     # plt.axhline(y=lower_bound, color='m', linestyle='--')
        #     # plt.axhline(y=upper_bound, color='m', linestyle='--')

        #     # plt.show()

        #     plt.clf()

        #     print("  - Q = {:}, ".format(Q))

        #     # If the power is within the confidence region, save the filtered data
        #     if (power_to_verify < upper_bound) and (power_to_verify > lower_bound):
        #         good_filtered_data = np.copy(filtered_data)

        #     # If the power is still outside the confidence region, increase the quality factor
        #     else:
        #         filtered_data = np.copy(good_filtered_data)
        #         Q += mains_hum[i]/config.mains_frequency
        # ### ###


        i += 1

    print("All frequencies filtered!")

    print()
    print("Finding background level")

    # Begin another least-squares optimization fit that follows the background level in the data
    # Begin initial guess's, 0 where additions/subtractions are and 1 where multiplications/divisions are
    x0 = np.array([1.0, 1.0])

    # Attempt to fit the background model on the filtered data
    res_lsq_back = scipy.optimize.least_squares(linearBackgroundModel, x0, 
                                                args=(time_relative, good_filtered_data)
                                                )

    print("Data background fit: {}".format(res_lsq_back.x))

    # Calculate the mean power using the least-squares solution
    background_level_arr = linearBackgroundModel(res_lsq_back.x, time_relative)


    ### Try to locate the fireball in the data and determine its boundaries ###

    if cml_args.integrate:

        # Calculate number of samples in window, since the main harmonice is 60/50 Hz which is 1/freq s, that will be the smallest window we should use 
        window_width = int(fs/mains_hum[0])
        window_shift = int(window_width/2.0)

        # To find the area beneath the curve, find the highest point in the data becuase it's likely the meteor
        peak_index = np.argmax(good_filtered_data)

        i = peak_index
        while(True):
            if(good_filtered_data[i]<background_level_arr[i] or i==0):
                left_std_start = i
                break
            else:
                i -= 1

        i = peak_index
        while(True):
            if(good_filtered_data[i]<background_level_arr[i] or i==len(good_filtered_data)-1):
                right_std_start = i
                break
            else:
                i += 1

        left_back_low = background_level_arr - 3*np.std(good_filtered_data[0:left_std_start])
        left_back_upp = background_level_arr + 3*np.std(good_filtered_data[0:left_std_start])

        right_back_low = background_level_arr - 3*np.std(good_filtered_data[right_std_start:-1])
        right_back_upp = background_level_arr + 3*np.std(good_filtered_data[right_std_start:-1])

        i = 0
        window_mean = np.mean(good_filtered_data[peak_index - (i + 1)*window_shift:peak_index - (i - 1)*window_shift])

        while(True):
            i += 1
            window_mean = np.mean(good_filtered_data[peak_index - (i + 1)*window_shift:peak_index - (i - 1)*window_shift])
            try:
                if(window_mean < left_back_upp[peak_index - i*window_shift] and window_mean > left_back_low[peak_index - i*window_shift]):
                    peak_lower_bound = peak_index - i*window_shift
                    break
            except IndexError:
                break

        i = 0
        window_mean = np.mean(good_filtered_data[peak_index + (i - 1)*window_shift:peak_index + (i + 1)*window_shift])

        while(True):
            i += 1
            window_mean = np.mean(good_filtered_data[peak_index + (i - 1)*window_shift:peak_index + (i + 1)*window_shift])
            try:
                if(window_mean < right_back_upp[peak_index + i*window_shift] and window_mean > right_back_low[peak_index + i*window_shift]):
                    peak_upper_bound = peak_index + i*window_shift
                    break
            except IndexError:
                break

        try:
            filtered_curve_area = scipy.integrate.trapz(good_filtered_data[peak_lower_bound:peak_upper_bound], time_relative[peak_lower_bound:peak_upper_bound])
            background_area = scipy.integrate.trapz(background_level_arr[peak_lower_bound:peak_upper_bound], time_relative[peak_lower_bound:peak_upper_bound])
            light_curve_relative_to_background = filtered_curve_area - background_area
        except NameError:
            print("Integration bounds could not be found.")


    ### ###

    

    ### Plot the spectrograms of the original data and the filtered data ###

    # Compute the spectrogram
    spectrum1, freqs1, t1, im1 = plt.specgram(intensity, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=0)
    spectrum2, freqs2, t2, im2 = plt.specgram(good_filtered_data, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=0)
    plt.close()


    min_val = 10 * np.log10(min(spectrum1.min(), spectrum2.min()))
    max_val = 10 * np.log10(min(spectrum1.max(), spectrum2.max()))

    fig, axs = plt.subplots(2, sharex=True, sharey=True)

    # Axis for top plot
    plt.sca(axs[0])
    axs[0].set_position([0.15, 0.5, 0.85, 0.4])
    axs[0].tick_params(bottom=False, labelbottom=False)
    # Axis for bottom plot
    plt.sca(axs[1])
    axs[1].set_position([0.15, 0.1, 0.85, 0.4])
    axs[1].tick_params(bottom=False, labelbottom=False)

    # Plotting original spectrogram on top plot
    plt.sca(axs[0])
    axs[0].set_title("Spectrogram before filtering")
    axs[0].set_ylabel("Frequency (Hz)")
    spectrum1, freqs1, t1, im1 = axs[0].specgram(intensity, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=0, vmin=min_val, vmax=max_val)


    # Plotting filtered spectrogram on bottom plot
    plt.sca(axs[1])
    axs[0].set_title("Spectrogram after filtering")
    axs[1].set_xlabel("Time")
    axs[1].set_ylabel("Frequency (Hz)")
    spectrum2, freqs2, t2, im2 = axs[1].specgram(good_filtered_data, Fs=sps, cmap='inferno', detrend='linear', NFFT=256, noverlap=0, vmin=min_val, vmax=max_val)
    
    cbar = fig.colorbar(im1, ax=axs.ravel().tolist(), pad=0.1, aspect = 30)
    cbar.ax.set_ylabel('Power (dB)')

    ax3 = axs[1].twiny()
    ax3.tick_params(bottom=True, labelbottom=True, top=False, labeltop=False)

    ax0_x0y0widthheight = list(axs[0].get_position().bounds) 
    ax1_x0y0widthheight = list(axs[1].get_position().bounds) 
    ax3_x0y0widthheight = list(ax3.get_position().bounds)

    ax3.set_position(ax1_x0y0widthheight)

    ax3.set_xlim(all_datetime[0],all_datetime[-1])

    plt.sca(ax3)
    plt.xticks(rotation=30)

    plt.savefig(os.path.join(csv_path, event_name + "_spectrograms.png"), dpi=150)

    if cml_args.showplots:
        plt.show()
    else:
        plt.clf()
        plt.close()


    

    
    # Plot the raw data, the filtered data, and the background level
    plt.plot(all_datetime, intensity, label='Raw Data')
    plt.plot(all_datetime, good_filtered_data, label='Filtered Data')
    plt.plot(all_datetime, background_level_arr, label='Background Level', color='green')

    if cml_args.integrate:

        plt.plot(all_datetime[0:left_std_start], left_back_upp[0:left_std_start], linestyle = '--', color = 'm')
        plt.plot(all_datetime[0:left_std_start], left_back_low[0:left_std_start], linestyle = '--', color = 'm')
        plt.plot(all_datetime[right_std_start:-1], right_back_upp[right_std_start:-1], linestyle = '--', color = 'm')
        plt.plot(all_datetime[right_std_start:-1], right_back_low[right_std_start:-1], linestyle = '--', color = 'm')
        try:
            plt.axvline(x=all_datetime[peak_upper_bound], color='red')
            plt.axvline(x=all_datetime[peak_lower_bound], color='red')
            plt.fill_between(all_datetime[peak_lower_bound:peak_upper_bound], good_filtered_data[peak_lower_bound:peak_upper_bound], background_level_arr[peak_lower_bound:peak_upper_bound], color='magenta')
        except NameError:
            pass

    plt.xticks(rotation=30)
    plt.legend()

    plt.savefig(os.path.join(csv_path, event_name + "_data.png"), dpi=150)

    if cml_args.showplots:
        plt.show()
    else:
        plt.clf()
        plt.close()

    corrected_datetime = all_datetime[int(0.02*len(all_datetime)):-1]
    corrected_good_data = good_filtered_data[int(0.02*len(good_filtered_data)):-1]
    corrected_background_level_arr = background_level_arr[int(0.02*len(background_level_arr)):-1]

    plt.plot(corrected_datetime, corrected_good_data, label='Filtered Data')
    plt.plot(corrected_datetime, corrected_background_level_arr, label='Background Level',color='green')
    plt.xticks(rotation=30)
    plt.title("Filtered Data")
    plt.legend()

    plt.savefig(os.path.join(csv_path, event_name + "_filtered_data.png"), dpi=150)

    if cml_args.showplots:
        plt.show()
    else:
        plt.clf()
        plt.close()

    plt.psd(intensity, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    plt.psd(good_filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    plt.title("Before and after filtering PSD")
    plt.plot(freqs,power_means, 'g')
    plt.plot(freqs,power_lower_lim, linestyle = '--', color = 'm')
    plt.plot(freqs,power_upper_lim, linestyle = '--', color = 'm')
    
    plt.savefig(os.path.join(csv_path, event_name + "_psd_filtered.png"), dpi=150)

    if cml_args.showplots:
        plt.show()
    else:
        plt.clf()
        plt.close()
