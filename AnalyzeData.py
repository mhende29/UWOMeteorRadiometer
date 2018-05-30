

from __future__ import print_function, division, absolute_import

import argparse
import os

import numpy as np
import scipy.signal
import scipy.optimize


DATA_SIZE = 2**20


class RDM(object):
    
    def __init__(self):
        
        self.header_size = None
        self.format_file_version = None
        self.station_code = None
        self.channel = None
        
        self.station_latitude = None 
        self.station_longitude = None
        self.station_elevation = None
        self.instrument_string = None
        
        self.num_samples = None
        self.checksum = None
        
        self.unix_start_s = None;
        self.unix_start_us = None;
        self.unix_end_s = None;
        self.unix_end_us = None;
        
        self.unix_s = None
        self.unix_us = None
        self.intensity = None


def readRDM(file_name, checksum_check=False, read_data=True):
    """ Read the binary RDM file. """

    with open(file_name, 'rb') as fid:
        
        rdm = RDM()
        
        rdm.header_size = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        rdm.format_file_version=np.fromfile(fid, dtype=np.uint32, count=1)[0]
        
        rdm.station_code = (b''.join(np.fromfile(fid, dtype='c', count=6))).decode("utf-8")
        # Reading 2 bytes even though there is only 1 char / there is an empty byte after the char
        rdm.channel = (b''.join(np.fromfile(fid, dtype='c', count=2))).decode("utf-8")
        
        rdm.station_latitude=np.fromfile(fid, dtype=np.double, count=1)[0]
        rdm.station_longitude=np.fromfile(fid, dtype=np.double, count=1)[0]
        rdm.station_elevation=np.fromfile(fid, dtype=np.double, count=1)[0]
        
        rdm.instrument_string = (b''.join(np.fromfile(fid, dtype='c', count=64))).decode("utf-8")
        rdm.num_samples = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        
        # Need to skip 4 bytes and not sure why
        np.fromfile(fid, dtype='c', count=4)
        
        rdm.checksum = np.fromfile(fid, dtype=np.uint64, count=1)[0]
        
        # Read file begin/end time
        rdm.unix_start_s = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        rdm.unix_start_us = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        rdm.unix_end_s = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        rdm.unix_end_us = np.fromfile(fid, dtype=np.uint32, count=1)[0]
        
        # Read the data arrays
        if read_data:
        
            # Skip to data
            fid.seek(rdm.header_size)
            
            # Read the tabular data
            table = np.fromfile(fid, dtype=np.uint32, count=3*(rdm.num_samples))
            table = np.reshape(table, (3, rdm.num_samples))
            
            # Unpack the values from the table
            rdm.time_s, rdm.time_us, rdm.intensity = table
            
            
            if checksum_check:
                
                # Compute the checksum
                chksum = np.sum(rdm.intensity.astype(np.uint64))
                
                # Compare the computed checksum to the checksum in the file
                if chksum == rdm.checksum:
                    return rdm, True
                    
                else:
                    return rdm, False
                    
            
        return rdm, False

        


def initNotchFilter(sps, freq, band, ripple, order, filter_type):
    """ Initializes a noth filter with the given parameters.

    Arguments:
        sps: [float] Samples per second.
        freq: [float] Middle frequency of the filter in Hz.
        band: [float] Width of the band in Hz.
        ripple: [foat] Ripple in the bandpass.
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



def filterLP(data, sps, mains_freq, lowpass=True, filter_order=3, additional=None):
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
    ripple = 1.0

    # Generate filter parameters for all harmonics
    filters_params = []
    for i in range(int((sps/2)/mains_freq)):

        # Compute the current harmonic frequency
        f_har = mains_freq*(i + 1)

        # Set proper filter band width, depending on the harmonic number (first ones are wider)
        if i == 0:
            band_har = 10.0
        elif i == 1:
            band_har = 7.5
        else:
            band_har = 5.0

        filters_params.append([sps, f_har, band_har, ripple, filter_order, filter_type])


    if additional is not None:
        for freq, band in additional:
            filters_params.append([sps, freq, band, ripple, filter_order, filter_type])


    filtered_data = np.copy(data)

    # Filter data
    for filt_param in filters_params:

        print(filt_param)

        # Init the filter
        b, a = initNotchFilter(*filt_param)

        # Filter the data
        filtered_data = scipy.signal.lfilter(b, a, filtered_data)



    if lowpass:

        ### Apply a lowpass filter which will filter out everything above the grid frequency ###

        # Init the lowpass filter
        Wn = mains_freq/(sps/2.0)
        b, a = scipy.signal.butter(3, Wn)

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




def fitSine(time_arr, signal_arr, guess_freq=None):
    """ Fit a sine to the given data. The sine will be fitted to the frequency with the highest amplitude. """

    if guess_freq is None:

        # Perform FFT, assume uniform spacing
        ff = np.fft.fftfreq(len(time_arr), (time_arr[1] - time_arr[0]))
        Fyy = np.abs(np.fft.fft(signal_arr))

        # Guess the frequency (excluding the zero frequency "peak", which is related to offset)
        guess_freq = np.abs(ff[np.argmax(Fyy[1:]) + 1])

    # Guess the amplitude and the offset
    guess_amp = np.std(signal_arr)*2.**0.5
    print(guess_amp)
    guess_offset = np.mean(signal_arr)

    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, phi, c):  
        return A * np.sin(w*t + phi) + c


    plt.plot(time_arr, signal_arr)
    plt.plot(time_arr, sinfunc(time_arr, *guess))
    plt.show()

    # Fit the sine
    popt, pcov = scipy.optimize.curve_fit(sinfunc, time_arr, signal_arr, p0=guess)

    A, w, phi, c = popt

    f = w/(2.*np.pi)

    fitfunc = lambda t: A * np.sin(w*t + phi) + c

    return {"amp": A, "omega": w, "phase": phi, "offset": c, "freq": f, "period": 1./f, "fitfunc": fitfunc, \
        "maxcov": np.max(pcov), "rawres": (guess,popt,pcov)}



if __name__ == "__main__":
        
    import matplotlib
    import matplotlib.pyplot as plt

    
    # Set up input arguments
    arg_p = argparse.ArgumentParser(description="Analyzes radiometer files.")
    arg_p.add_argument('rdm_file', type=str, help="Path to the .rdm file.")
    
    # Parse input arguments
    cml_args = arg_p.parse_args()
    
    
    # Read the binary RDM file
    rdm, chksum_pass = readRDM(cml_args.rdm_file)
    
    # Tell us if the chksum passed
    print(chksum_pass)
    
    # Print header data
    print(rdm.header_size)
    print(rdm.format_file_version)
    print(rdm.station_code)
    print(rdm.channel)
    print(rdm.station_latitude)
    print(rdm.station_longitude)
    print(rdm.station_elevation)
    print(rdm.instrument_string)
    print(rdm.num_samples)
    print(rdm.checksum)
    
    print(rdm.unix_start_s)
    print(rdm.unix_start_us)
    print(rdm.unix_end_s)
    print(rdm.unix_end_us)
        
    # Print the tabular data
    print(rdm.intensity)
    print(rdm.time_s)
    print(rdm.time_us)
    
    # Convert UNIX times from int to one float
    unix_times = rdm.time_s.astype(np.float64) + rdm.time_us.astype(np.float64)/1e6
    print(unix_times)

    # Compute relative time since the beginning of the recording
    time_relative = unix_times - np.min(unix_times)
    
    # Compute samples per second
    sps = len(rdm.time_s)/(unix_times[-1] - unix_times[0])

    print('SPS:', sps)
    
    # Extract file name
    _, file_name = os.path.split(cml_args.rdm_file)
    file_name = file_name.replace('.rdm', '')


    # # Plot raw data
    # plt.plot(time_relative, rdm.intensity)
    # plt.show()
    

    # Plot the spectrogram
    plt.specgram(rdm.intensity, Fs=sps, cmap='inferno', detrend='linear', NFFT=2048, noverlap=0)
    
    plt.title(file_name)
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')
    
    # plt.savefig('/home/pi/Desktop/{:s}.png'.format(file_name), dpi=300)
    
    plt.show()


    # Plot power spectral density
    plt.psd(rdm.intensity, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    plt.show()


    # Filter the light pollution
    mains_freq = 60.0 # Hz
    filtered_data = filterLP(rdm.intensity, sps, mains_freq, additional=[(20, 2.0), (32, 2.0), (94, 2.0), (553.0, 2.0), (614, 2.0)], lowpass=False)


    # # Apply a moving average
    # window_size = 3
    # filtered_data = movingAverage(rdm.intensity, n=window_size)
    # time_relative = time_relative[:len(filtered_data)]
    # filtered_data = filtered_data[::window_size]
    # time_relative = time_relative[::window_size]
    # sps = sps/window_size


    # # Apply a broad lowpass and a highpass filter
    # bandpass_low = 5.0 # Hz
    # bandpass_high = mains_freq
    # filtered_data = filterBandpass(filtered_data, sps, bandpass_low, bandpass_high, order=3)

    print(filtered_data)

    # # Fit a sine
    # sine_fit = fitSine(time_relative[int(len(time_relative)*0.1):], filtered_data[int(len(time_relative)*0.1):], guess_freq=0.5)

    # plt.plot(time_relative, filtered_data)
    # plt.plot(time_relative, sine_fit['fitfunc'](time_relative))

    # plt.show()

    # #filtered_data -= sine_fit['fitfunc'](time_relative)


    # Plot filtered spectrogram
    plt.specgram(filtered_data, Fs=sps, cmap='inferno', detrend='linear', NFFT=2048, noverlap=0)

    plt.title(file_name + ' filtered')
    plt.xlabel('Time (s)')
    plt.ylabel('Frequency (Hz)')

    plt.show()


    # Plot power spectral density
    plt.psd(filtered_data, Fs=sps, detrend='linear', NFFT=2048, noverlap=0)
    plt.show()


    # Plot filtered data
    plt.plot(time_relative, filtered_data)
    plt.show()
    
    

    
    
    
    
