
from __future__ import print_function

import os
import sys
import argparse
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime

from getRDMData import readRDM
from AnalyzeData import movingAverage


def getNightPlot(file_path):
	""" Generates two plots of the nights data. 

		Input Arguments:
			
			-file_path (string): The path to the directory that stores the nights data. Ex: /home/user/Desktop/RadiometerData/CaptureData/AA0000_A_19700101-000000

		Outputs:

			- Two .png plots saved in the given file_path
	"""

	# Enter the nights directory and gather the names of all the .rdm files 
	file_list = [file_name for file_name in os.listdir(file_path) if file_name.endswith(".rdm")]
	file_list = sorted(file_list)

	# Define how many .rdm files were in the directory
	file_size = len(file_list)
	files_left = file_size

	# Define several empty lists to store the data that will be collected
	data_average = []
	data_peaks = []
	data_times = []
	raw_data = []
	raw_time = []
	zeros_index = []

	# We have 4 MB of data per file which translates to 1048576 or 2^20, since we're taking bins 64 samples which is 2^6 samples
	# we will need to run 2^20/2^6 loops which is 2^(20-6) = 2^14 loops per file
	bin_size = 8192 
	loops = int(2**20/bin_size)
	
	# Gather data from the files one at a time to reduce RAM usage
	for file_name in file_list:
		
		# Read the RDM file
		rdm, rdm_status = readRDM(os.path.join(file_path, file_name))

		# Combine the micro-seconds and the seconds
		time = rdm.time_s.astype(np.float64) + rdm.time_us.astype(np.float64)/1e6

		# Compute the moving average
		intensity = movingAverage(rdm.intensity, 256)



		# Compute the data averages and peaks by binning in order to reduce the amount of data plotted by 8192 times
		for i in range(0,loops):

			temp_data = intensity[i*bin_size:(i+1)*bin_size]
			temp_time = time[i*bin_size:(i+1)*bin_size]

			if(len(temp_data) != 0):
				data_peaks.append(temp_data.max())
				data_average.append(np.average(temp_data))
				data_times.append(np.average(temp_time))
			
		files_left -= 1

		# Update on the averaging process
		print("{:.2%} Complete".format(1.0-files_left/float(file_size)))

	# Convert the averages and peaks to arrays	
	data_average = np.array(data_average)
	data_peaks = np.array(data_peaks)

	# Create a list of datetime objcets from the unix times
	data_times = [datetime.utcfromtimestamp(t) for t in data_times]

	# Plot the averages and the peaks
	plt.plot(data_times, data_average, label = "Averages")
	plt.plot(data_times, data_peaks, label = "Peaks")
	plt.grid(which = "both")
	plt.gca().set_yscale('log')
	plt.xticks(rotation=30)
	plt.xlabel('Time')
	plt.ylabel('ADU')
	plt.legend(loc = 0)

	# Determine the appropriate title
	if(data_times[0].strftime('%d') == data_times[-1].strftime('%d')):
		plt.title('Peak intensities for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d, %Y'))
	else:
		plt.title('Peak intensities for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d') + ' to ' + data_times[-1].strftime('%d, %Y'))

	# Save the figure
	plt.savefig(os.path.join(file_path, "NightPlot_{:s}.png".format(data_times[0].strftime('%Y%m%d'))), dpi=300)
	plt.close()

	# Determine if 
	zeros_index = [i for i in range(len(data_average)) if(data_average[i] == 0)]
	
	# Plot the max/minus of the data (shows where peaks are more clearly)
	max_minus = data_peaks[1:-1] - (data_average[:-2]+data_average[2:])/2

	# If for some reason there is garbage data in the average set the max_minus at those locations to 0, +/-2 to cover the sides as well
	if (len(zeros_index)!=0):
		for zeros in zeros_index:
			max_minus[zeros-2:zeros+2] = 0

	plt.plot(data_times[1:-1], max_minus)
	plt.grid(which = "both")
	plt.xticks(rotation=30)
	plt.xlabel('Time')
	plt.ylabel('ADU')
	
	# Determine the appropriate title
	if(data_times[0].strftime('%d') == data_times[-1].strftime('%d')):
		plt.title('Intensity variations for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d, %Y'))
	else:
		plt.title('Intensity variations for ' + file_list[0][0:6] + ' on the night of ' + data_times[0].strftime('%B %d') + ' to ' + data_times[-1].strftime('%d, %Y'))
	
	# Save the figure
	plt.savefig(os.path.join(file_path, "MaxMinus_{:s}.png".format(data_times[0].strftime('%Y%m%d'))), dpi=300)
	plt.clf()
	
if __name__ == "__main__":
	
	# Set up input arguments
	arg_p = argparse.ArgumentParser(description="Get radiometer files.")
	arg_p.add_argument('dir_path', type=str, help="Path to the .rdm file.")

	# Parse input arguments
	cml_args = arg_p.parse_args()
	
	getNightPlot(cml_args.dir_path)
	print("Done")
	
