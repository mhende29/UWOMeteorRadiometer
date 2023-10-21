
""" Gathers all rdm files within a certain directory, and searchs for the two files that include data near the given time. 
	
	Input arguments: 
						-file_path: The """


import os
import sys
import argparse
import tarfile
import shutil
import numpy as np
from datetime import datetime, timedelta

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
	""" Read the binary RDM file. 

		Input Arguments:
			
			-file_name (string): The path and name of the file we want to decode. Ex: /home/user/Desktop/RadiometerData/ArchivedData/AA0000_A_19700101-000000.000000_000000.000000.rdm
			-checksum_check (boolean): If enabled compares the sum of all intensities to the value stored before compression to ensure data is valid
			-read_data (boolean): If enabled outputs all data, if not outputs file information

		Outputs:

			- rdm (rdm object): The rdm object containing all data stored in the rdm file
			- rdm_status (boolean): Boolean indicating wether the data is valid or not, ignore if checksum_check = False
	"""

	with open(file_name, 'rb') as fid:
		
		# Define an rdm object
		rdm = RDM()
		
		# Load the header size and the formatting version
		rdm.header_size = np.fromfile(fid, dtype=np.uint32, count=1)[0]
		rdm.format_file_version=np.fromfile(fid, dtype=np.uint32, count=1)[0]

		# Read the station code and the station channel
		rdm.station_code = (b''.join(np.fromfile(fid, dtype='c', count=6))).decode("utf-8")
		# Reading 2 bytes even though there is only 1 char / there is an empty byte after the char
		rdm.channel = (b''.join(np.fromfile(fid, dtype='c', count=2))).decode("utf-8")

		# Load the geographical coordinates
		rdm.station_latitude=np.fromfile(fid, dtype=np.double, count=1)[0]
		rdm.station_longitude=np.fromfile(fid, dtype=np.double, count=1)[0]
		rdm.station_elevation=np.fromfile(fid, dtype=np.double, count=1)[0]

		# Load the description of the device
		rdm.instrument_string = (b''.join(np.fromfile(fid, dtype='c', count=64))).decode("utf-8")
		rdm.num_samples = np.fromfile(fid, dtype=np.uint32, count=1)[0]
		
		# Need to skip 4 bytes and not sure why
		np.fromfile(fid, dtype='c', count=4)
		
		# Load the checksum
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


			# Skip all values where the UNIX time (seconds) is 0
			rdm.time_s, rdm.time_us, rdm.intensity = table[:, rdm.time_s > 0]
					
		return rdm, False

        
def unzipData(night_path, desired_files):
	""" Unzips all desired files and saves their data in two lists. 

		Input Arguments:
			
			-night_path (string): The path to the list of files we want to decode. Ex: /home/user/Desktop/RadiometerData/ArchivedData
			-desired_files (list of strings): List containing the names of all the zipped files we want to read. Ex: AA0000_A_19700101-000000.000000_000000.000000.rdm.tar.bz2

		Outputs:

			- all_data (list of ints): The rdm intensities in the desired files
			- all_time (list of floats): The rdm timings in the desired files
	"""

	# Create an empty list to store the rdm objects
	rdm_obj = []

	# Get working directory
	temp_dir = os.path.join(os.getcwd(),"temp")

	# If there is no temporary working directory make one
	if(os.path.isdir(temp_dir)):
		pass
	else:
		os.mkdir(temp_dir ,0o755)

	# Unzip all desired files, one at a time
	for zipped_file in desired_files:

		# If the file is not zipped, load it directly
		if zipped_file.endswith('.rdm'):

			# Read the RDM file, and add the object to the list
			temp_rdm, rdm_status = readRDM(os.path.join(night_path, zipped_file))
			rdm_obj.append(temp_rdm)

			continue



		with tarfile.open(os.path.join(night_path,zipped_file),"r:bz2") as tar_file:

			# Get the member and name of each zipped file we go through, getmembers and getnames return a list of 1 
			member = tar_file.getmembers()
			name = tar_file.getnames()

			# Extract the current desired file
			tar_file.extractall(temp_dir,member)

			# Read the RDM file, and add the object to the list
			temp_rdm, rdm_status = readRDM(os.path.join(temp_dir, *name))
			rdm_obj.append(temp_rdm)

			# Delete the unpacked file once the data has been gathered
			os.remove(os.path.join(temp_dir, *name))

	# Remove the temporary working directory
	shutil.rmtree(temp_dir)
	
	# Create a list to store all time and intensity data
	all_data, all_time = [], []

	# Append the data into two arrays containing all required information
	for i in range(len(rdm_obj)):
		all_data.append(rdm_obj[i].intensity.astype(np.int))
		all_time.append(rdm_obj[i].time_s.astype(np.float64) + rdm_obj[i].time_us.astype(np.float64)/1e6)

	return all_data, all_time




def getRDMData(dir_path, station_code, station_channel, time, time_range, UT_corr = 0.0):

	""" Determines which data needs to be read and retrieves the data. 

		Input Arguments:
			
			-dir_path (string): The path to the directory storing all zipped files. Ex: /home/user/Desktop/RadiometerData/ArchivedData
			-station_code (string): The station code of the data we want to read. Ex: AA0000
			-station_channel (char): The station channel of the data we want to read. Ex: A
			-time (string): The time period at the centre of the analysis, follows this format: YYYYMMDD-HHmmss.uuuuuu & YYYYMMDD-HHmmss Ex: 20180703-040000.06 or 20180703-040000
			-time_range (float): Number of seconds the total analysis covers with "time" being at the centre of the time range. Ex: 500
			-UT_corr (float): UTC correction time in hours

		Outputs:

			- valid_data (list of ints): The rdm intensities in the desired time section
			- valid_time (list of floats): The rdm timings in the desired time section
	"""	

	# Create a datetime object from the given time
	# First try if the given time ends in microseconds
	try:
		desired_time = datetime.strptime(time,"%Y%m%d-%H%M%S.%f") - timedelta(hours = UT_corr)
	
	# If there is a value error because there wasn't microseconds try without it
	except ValueError:
		desired_time = datetime.strptime(time,"%Y%m%d-%H%M%S") - timedelta(hours = UT_corr)
		
	# Adds and subtracts half the time range in seconds to give bounds to the data
	begin_time = desired_time - timedelta(0,time_range/2)
	end_time = desired_time + timedelta(0,time_range/2)

	# Get station code and channel in order to find valid zipped files
	station = station_code + "_" + station_channel + "_"
	
	# Gather the directories of all the nights that contain the station code and channel
	potential_nights = [night_name for night_name in os.listdir(dir_path) if night_name.startswith(station)]

	# Sorts the files alphabetically which also sorts them by age
	potential_nights = sorted(potential_nights)

	# Creates a list with an identical length to the potential zip folder however it store the startime of each zip file as a datetime object
	night_start_times = [datetime.strptime(night_name,station + "%Y%m%d-%H%M%S") - timedelta(hours = UT_corr) for night_name in potential_nights]

	# Check to see if the user asked to see data older than the earliest data.
	if(night_start_times[0] > desired_time):
		print("Earliest files on record began:", night_start_times[0], ".\nPlease provide a valid time period that occured after the oldest period.")
		sys.exit()

	# Create an index to find the desired file
	index = 0

	# Search for the datetime where our desired datetime is found between two zips
	for time in night_start_times:

		if(time == night_start_times[0]):
			prev_time = time

		# If the date were looking at is after the desired date and the previous one was before the desired date, we found the desired date
		if(time > desired_time and prev_time < desired_time):
			# Correct the index since the file we want is associated with the previous time
			index -= 1
			break
		
		# Desired file not found yet, increment index
		index += 1

		# Keep track of the previous datetime
		prev_time = time


	# Still havent yet found a file that contains the data we want, can't be guaranteed the file exists because we only have a lower bound to compare to
	# Therefore check the final folder to be safe
	if(index == len(potential_nights)):
		index -= 1

	# Found the directory for the desired nights data
	desired_dir = potential_nights[index]

	# Gets the path to the desired nights data
	night_path = os.path.join(dir_path,desired_dir)

	# Creates a list of the zipped files in the nights directory and sorts them by age
	rdm_list = [file_name for file_name in os.listdir(night_path) if (file_name.endswith(".tar.bz2") or file_name.endswith(".rdm"))]
	rdm_list = sorted(rdm_list)

	# Get the start and end times of each file
	rdm_start_times = [datetime.strptime(file_name, station + "%Y%m%d-%H%M%S.%f" + file_name[31:]) - timedelta(hours = UT_corr) for file_name in rdm_list]
	rdm_end_times = [datetime.strptime(file_name, station + "%Y%m%d" + file_name[17:32] + "%H%M%S.%f" + file_name[45:]) - timedelta(hours = UT_corr) for file_name in rdm_list]
	
	# Check to see if the user asked to see data older than the earliest data.
	if(rdm_end_times[-1] < desired_time):
		print("Latest file on record stopped at:",rdm_end_times[-1],".\nPlease provide a valid time period that occured before the most recent period.")
		sys.exit()

	# Create an empty list to store the desired files in the nights directory
	desired_files = []

	# If the files contain data within the defined bounds, remember the name of the files
	for i in range(len(rdm_list)):
		if((rdm_end_times[i] > begin_time) and (rdm_start_times[i] < end_time)):
			desired_files.append(rdm_list[i])

	all_data, all_time = unzipData(night_path, desired_files)

	# Convert the nested lists into numpy arrays
	all_data = np.array([data for data_file in all_data for data in data_file])
	all_time = np.array([time for data_file in all_time for time in data_file])

	# Create a list containing all times as datetime objects
	all_datetime = [datetime.utcfromtimestamp(t) for t in all_time]

	# Get the unix equivalent of the time bounds
	desired_lower_unix = (begin_time - datetime(1970, 1, 1)).total_seconds()
	desired_upper_unix = (end_time - datetime(1970, 1, 1)).total_seconds()

	# Find indices to cut
	beg_ind = (np.abs(all_time - desired_lower_unix)).argmin()
	end_ind = (np.abs(all_time - desired_upper_unix)).argmin()

	# Store the valid data into three arrays
	valid_data = all_data[beg_ind:end_ind]
	valid_time = all_time[beg_ind:end_ind]
	valid_datetime = all_datetime[beg_ind:end_ind]

	return valid_data, valid_time


if __name__ == "__main__":
	
	# Set up input arguments
	arg_p = argparse.ArgumentParser(description="Get radiometer files.")
	arg_p.add_argument('time', metavar='TIME', nargs='?', \
		help='Time of the event in the YYYYMMDD-HHMMSS.ms format.', type=str, default=None)
	arg_p.add_argument('range', metavar='DURATION_SECONDS', help="""Grabs data so the total, 
		range of data covers the range. Range/2 on both sides of time. """)
	
	# Parse input arguments
	cml_args = arg_p.parse_args()

	file_stored = "/home/michael/Desktop/Radiometer Events/ArchiveTest"
	station_code = "CA0001"
	station_channel = "A"

	getRDMData(file_stored, station_code, station_channel, cml_args.time, cml_args.range)
	print("Done")