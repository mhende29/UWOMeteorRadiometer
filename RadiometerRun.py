

from __future__ import print_function, absolute_import

import os
import time
import math
import datetime
import argparse
import configparser
import ads1256
import signal
import sys
import multiprocessing

from threading import Event
from CaptureDuration import captureDuration
from DeleteOldObservations import deleteOldObservations
from Misc import archiveDir
from UploadManager import UploadManager
from NightPlot import getNightPlot
from GetRDMConfig import RDMConfig, readConfig, makeConfig

# Flag indicating that capturing should be stopped
STOP_CAPTURE = False

def breakHandler(signum, frame):
    """ Handles what happens when Ctrl+C is pressed. """
    
    print("\nProgram forced stopped via ctrl+c!\n")
    global STOP_CAPTURE

    # Set the flag to stop capturing video
    STOP_CAPTURE = True
    
    exit_loop.set()
    exit_loop.clear()
    
    

# Save the original event for the Ctrl+C
ORIGINAL_BREAK_HANDLE = signal.getsignal(signal.SIGINT)


def setSIGINT():
    """ Set the breakHandler function for the SIGINT signal, will be called when Ctrl+C is pressed. """

    signal.signal(signal.SIGINT, breakHandler)

def resetSIGINT():
    """ Restore the original Ctrl+C action. """

    signal.signal(signal.SIGINT, ORIGINAL_BREAK_HANDLE)


def timeToRun(time_hrs):
    run_seconds = 500*math.ceil(7.2*time_hrs)
    hours = math.floor(run_seconds/3600)
    minutes = math.floor((run_seconds/3600-hours)*60)
    seconds = math.floor(((run_seconds/3600-hours)*60 - minutes)*60)
    return hours,minutes,seconds

def makeDirectory(station_code,channel,path):
    
    # Get the current time
    creation_time = datetime.datetime.utcnow()
    
    # Generate the name of the directory based on configuration information and the time of excution
    directory_name = station_code + "_" + channel + '_{:04d}{:02d}{:02d}-{:02d}{:02d}{:02d}'.format(creation_time.year,creation_time.month,creation_time.day,creation_time.hour,creation_time.minute,creation_time.second)

    # Create the full path name
    directory = os.path.join(path, directory_name)
    
    # Create the directory to store current data 
    os.mkdir(directory,0o755) 
    
    # Returns the path where the c code can store the data, also adds a slash at the tail of the directory as
    # the c code requires it and double the amount of slashes as c also requires it.
    return (directory + "/").replace("/","//")
    
class UploadConfig(object):
    def __init__(self,stationID,hostname,remote_dir,rsa_private_key,upload_queue_file):

        # remote hostname where SSH server is running
        self.stationID = "rdm_" + stationID.lower()         
        self.hostname = hostname
        self.remote_dir = remote_dir
        self.rsa_private_key = os.path.expanduser(rsa_private_key)
        self.upload_queue_file = upload_queue_file
        
##############################################################################################################

exit_loop = Event()

if __name__ == "__main__":
    
    # Since this code is to run on a raspberry pi, define a default home for all data, and a default name
    # for the configuration file.
    DEFAULT_PATH = "/home/pi/RadiometerData"
    CONFIG = "config.txt"
    
    setSIGINT()
    
    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description=""" Starting capture and compression.""")
    
    # Add a mutually exclusive for the parser (the arguments in the group can't be given at the same time)
    arg_group = arg_parser.add_mutually_exclusive_group()
    
    # Duration and loops are both mutually exclusive because they both provide a length of time to run for
    arg_group.add_argument('-d', '--duration', metavar='DURATION_HOURS', help="""Start capturing right away, 
        with the given duration in hours. """)
    arg_group.add_argument('-l', '--loops', metavar='NUMBER_LOOPS', help="""Start capturing right away, 
        with the given duration in loops. """)
        
    # check to see if a config file was provided
    arg_parser.add_argument('-c', '--config', metavar="FILE_PATH", help="""Read the given config file as, 
        opposed to read the file from the default path. """)
    
    # Parse the command line arguments
    cml_args = arg_parser.parse_args()
    
    # Check if the home directory exists, if not create it and exit
    if (not os.path.isdir(DEFAULT_PATH)):
        os.mkdir(DEFAULT_PATH,0o755)
        print("A new directory called RadiometerData has been created to store all data collected.\nIf no config.txt is found in this directory on the next run,\none will be generated and can then be generated and can be altered.\n\nThe program will now close.")
        sys.exit()
    
    # Check to see if a default config file was given
    if cml_args.config:
        rdm_config = readConfig(cml_args.config)
        
    # Check to see if there is a config file
    elif (os.path.isfile(os.path.join(DEFAULT_PATH,CONFIG))):
        rdm_config = readConfig(os.path.join(DEFAULT_PATH,CONFIG))
        
    else:
        makeConfig(os.path.join(DEFAULT_PATH,CONFIG))
    
    # Check if storage files exist, if not create them
    # Check if both don't exist
    if (not (os.path.isdir(os.path.join(DEFAULT_PATH, rdm_config.raw)) or os.path.isdir(os.path.join(DEFAULT_PATH, rdm_config.zipped)))):
        os.mkdir(os.path.join(DEFAULT_PATH, rdm_config.raw),0o755)
        os.mkdir(os.path.join(DEFAULT_PATH, rdm_config.zipped),0o755)
        
    # Check if the raw data folder exists and the archiving folder doesn't yet
    elif (os.path.isdir(os.path.join(DEFAULT_PATH, rdm_config.raw)) and (not os.path.isdir(os.path.join(DEFAULT_PATH, rdm_config.zipped)))):
        os.mkdir(os.path.join(DEFAULT_PATH, rdm_config.zipped),0o755)
        
    # Check if the archiving folder exists and the raw data folder doesn't yet
    elif (os.path.isdir(os.path.join(DEFAULT_PATH, rdm_config.zipped)) and (not os.path.isdir(os.path.join(DEFAULT_PATH, rdm_config.raw)))):
        os.mkdir(os.path.join(DEFAULT_PATH, rdm_config.raw),0o755)
    
    upload_config = UploadConfig(rdm_config.station_code, rdm_config.hostname, rdm_config.remote_dir, rdm_config.rsa_private_key, rdm_config.upload_queue_file)
    
        
##############################################################################################################
    
    # If the duration of capture was given, capture right away for a specified time in hours
    if cml_args.duration:
        
        # Convert the given duration to floating point
        duration = float(cml_args.duration)
    
    elif cml_args.loops:
        
        # Convert the given number of loops to floating point hours by dividing loops by loops per hour
        duration = (float(cml_args.loops))/7.2

    upload_manager = None
    if rdm_config.upload_enabled:
        # Init the upload manager
        upload_manager = UploadManager(upload_config)
        upload_manager.start()


##############################################################################################################

    # If a duration was given in loops or hours, run for the desired length of time
    if (cml_args.duration or cml_args.loops):
        
        deleteOldObservations(DEFAULT_PATH, rdm_config.raw, rdm_config.zipped, duration)
        
        stored_path = makeDirectory(rdm_config.station_code, rdm_config.channel, os.path.join(DEFAULT_PATH, rdm_config.raw))

        resetSIGINT()

        # ads1256.run returns 0 when its done recording
        running = ads1256.run(duration, rdm_config.mode, rdm_config.gain, rdm_config.station_code, rdm_config.channel, rdm_config.latitude, rdm_config.longitude, rdm_config.elevation, rdm_config.instrument_string, stored_path)

        setSIGINT()
        
        # Compress all files saved
        source = stored_path.replace("//","/")[:-1]
        
        # If the ADC had a faulty exit, end the program
        if running == 1:
            print("\nProgram terminated!\n")
        # If it exited safely, zip the data
        else:
            
            # Zips the files and moves the =m to the archiving directory
            archive_name = archiveDir(source, os.listdir(source), os.path.join(os.path.join(DEFAULT_PATH, rdm_config.zipped),os.path.split(source)[1]) , os.path.split(source)[1])
            
            # Put the archive up for upload
            if upload_manager is not None:
                upload_manager.addFiles([archive_name])
    
    
    
    
##############################################################################################################

    # If no duration given, it is assumed to wait for night and run all night.
    else:
        while(not STOP_CAPTURE):
            
            # Get a datetime object of exactly when sunset is at -5.5 deg and how long to run for until sunrise in seconds
            start_time, duration = captureDuration(rdm_config.latitude, rdm_config.longitude, rdm_config.elevation)
            
            # Convert duration to hours
            duration/=3600
            
            deleteOldObservations(DEFAULT_PATH, rdm_config.raw, rdm_config.zipped, duration)
            
            # If start_time holds a datetime object as opposed to being True
            if start_time != True:
                
                # Get a current datetime object
                time_now = datetime.datetime.utcnow()
                # Calculate the difference in time between now and when the recording should start
                waiting_time = start_time - time_now

                hours, minutes, seconds = timeToRun(duration)

                # Let the user know when it will start recording in UTC
                print("Recording will start on:", start_time, "UTC and\nrecord for {:02d} hours, {:02d} minutes and {:02d} seconds.".format(hours, minutes, seconds))
                
                # Wait until sunset using the wait time in seconds
                exit_loop.wait(int(waiting_time.total_seconds()))
                if (STOP_CAPTURE):
                    break
                
            # Creates the directory to store the data from this current run and returns the file path used to 
            # store the data in c, uses two slashes instead of one "/" -> "//"
            stored_path = makeDirectory(rdm_config.station_code, rdm_config.channel, os.path.join(DEFAULT_PATH, rdm_config.raw))
            
            resetSIGINT()
            
            # ads1256.run returns 0 when its done recording, or 1 for a forced exit, automatically saves data
            # to stored_path
            running = ads1256.run(duration, rdm_config.mode, rdm_config.gain, rdm_config.station_code, rdm_config.channel, rdm_config.latitude, rdm_config.longitude, rdm_config.elevation, rdm_config.instrument_string, stored_path)
            
            setSIGINT()
            
            # Remove the extra slash used at the end of the path (needs to be there for c code), and convert 
            # back to single slashes. The last slash is ignored by taking all but the last element
            source = stored_path.replace("//","/")[:-1]
            
            # If the ADC had a faulty exit, end the program
            if running == 1:
                print("\nProgram terminated with ctrl c!\n")
                break
                
            # If it exited safely, zip the data
            else:
                
                getNightPlot(source)
                
                # Zips the files and moves the =m to the archiving directory
                archive_name = archiveDir(source, os.listdir(source), os.path.join(os.path.join(DEFAULT_PATH, rdm_config.zipped),os.path.split(source)[1]) , os.path.split(source)[1])
                
                # Put the archive up for upload
                if upload_manager is not None:
                    upload_manager.addFiles([archive_name])
            
    # Stop the upload manager
    if upload_manager.is_alive():
        upload_manager.stop()
        del upload_manager
        
    sys.exit()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
