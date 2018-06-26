
from __future__ import print_function

import os
import sys
import argparse
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime

from AnalyzeData import readRDM, movingAverage


def getNightPlot(file_path, get_peaks = True):



	file_list = [file_name for file_name in os.listdir(file_path)]
	file_list = sorted(file_list)
	file_size = len(file_list)
	files_left = file_size

	data_average = []
	data_peaks = []
	data_times = []
	raw_data = []
	raw_time = []
	# We have 4 MB of data per file which translates to 1048576 or 2^20, since we're taking bins 64 samples which is 2^6 samples
	# we will need to run 2^20/2^6 loops which is 2^(20-6) = 2^14 loops per file
	bin_size = 8192 
	loops = int(2**20/bin_size)
	

	for file_name in file_list:
		
		# Read the RDM file
		rdm, rdm_status = readRDM(os.path.join(file_path, file_name))

		time = rdm.time_s.astype(np.float64) + rdm.time_us.astype(np.float64)/1e6

		# Compute the moving average
		intensity = movingAverage(rdm.intensity, 256)
		
		# plt.plot(time[63:], intensity)
		# plt.show()

		for i in range(0,loops):
			temp_data = intensity[i*bin_size:(i+1)*bin_size]
			temp_time = time[i*bin_size:(i+1)*bin_size]
			data_peaks.append(temp_data.max())
			data_average.append(np.average(temp_data))
			data_times.append(np.average(temp_time))
		files_left -= 1

		print("{:.2%} Complete".format(1.0-files_left/float(file_size)))
		
	#plt.plot(data_times,data_average)
	#plt.plot(data_times,data_peaks)
	data_average = np.array(data_average)
	data_peaks = np.array(data_peaks)
	data_times = [datetime.utcfromtimestamp(t) for t in data_times]

	plt.plot(data_times,data_average)
	plt.plot(data_times,data_peaks)
	plt.grid(which = "both")
	plt.gca().set_yscale('log')
	plt.xticks(rotation=30)
	plt.xlabel('Time')
	plt.ylabel('ADU')

	if(data_times[0].strftime('%d') == data_times[-1].strftime('%d')):
		plt.title('Peak intensities for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d, %Y'))
	else:
		plt.title('Peak intensities for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d') + ' to ' + data_times[-1].strftime('%d, %Y'))
	plt.show()

	plt.plot(data_times[1:-1], data_peaks[1:-1] - (data_average[:-2]+data_average[2:])/2)
	plt.grid(which = "both")
	plt.xticks(rotation=30)
	plt.xlabel('Time')
	plt.ylabel('ADU')
	if(data_times[0].strftime('%d') == data_times[-1].strftime('%d')):
		plt.title('Intensity variations for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d, %Y'))
	else:
		plt.title('Intensity variations for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d') + ' to ' + data_times[-1].strftime('%d, %Y'))
	plt.show()

if __name__ == "__main__":

	#file_stored = "/home/michael/Desktop/Radiometer Events/June 20 2018/"
	
	# Set up input arguments
	arg_p = argparse.ArgumentParser(description="Get radiometer files.")
	arg_p.add_argument('dir_path', type=str, help="Path to the .rdm file.")

	# Parse input arguments
	cml_args = arg_p.parse_args()
	
	getNightPlot(cml_args.dir_path)
	print("Done")
	