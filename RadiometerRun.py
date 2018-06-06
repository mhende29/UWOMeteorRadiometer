

from __future__ import print_function, absolute_import

import os
import time
import datetime
import argparse
import configparser
import ads1256
import signal
import sys

from collections import OrderedDict
from CaptureDuration import captureDuration
from Misc import archiveDir

def makeDirectory(station_code,channel,path):
    path = "/home/pi/RadiometerData/CapturedData/"
    
    creation_time = datetime.datetime.utcnow()
    
    directory_name = station_code + "_" + channel + '_{:04d}{:02d}{:02d}-{:02d}{:02d}{:02d}'.format(creation_time.year,creation_time.month,creation_time.day,creation_time.hour,creation_time.minute,creation_time.second)

    directory = path + directory_name
    
    os.mkdir(directory,0o755) 
    
    return (directory + "/").replace("/","//")

##############################################################################################################

if __name__ == "__main__":
    
    # Check to see if there is a config file
    if (os.path.isfile('/home/pi/RadiometerData/config.txt')):
        
        # Create a config object
        config = configparser.ConfigParser()
    
        # Read the config file into the object
        config.read('/home/pi/RadiometerData/config.txt')
    
        # Gather configuration data
        station_code = config['Station']['StationCode']
        channel = config['Station']['Channel']
        latitude = float(config['Station']['Latitude'])
        longitude = float(config['Station']['Longitude'])
        elevation = float(config['Station']['Elevation'])
        instrument_string = config['Station']['InstrumentString']
        path = config['Station']['Path']
        raw = config['Station']['RawData']
        zipped = config['Station']['StoredData']
        mode = int(config['Station']['DifferentialMode'])
        
    else:
        
        # There was no detected config file so one will be created
        # An error message explaining the issue
        print("No config file detected in /home/pi/RadiometerData/")
        print("A default config file has been created and can be changed in RadiometerData")
        
        # Create a config object
        config = configparser.ConfigParser()
        
        # optionxform prevents it from naming all config parameters with lower case letters
        config.optionxform = str
        
        # Creates the data inside the config file using default values
        config['Station'] = OrderedDict((
            ('StationCode', 'AA0000'),
            ('Channel', 'A'),
            ('Latitude', '0.0'),
            ('Longitude', '0.0'),
            ('Elevation', '0.0'),
            ('InstrumentString', 'Your description'),
            ('Path','/home/pi/RadiometerData/'),
            ('RawData','CapturedData'),
            ('StoredData','ArchivedData'),
            ('DifferentialMode','1')
        ))
        
        # Generate the file in the desired directory and close it
        with open('/home/pi/RadiometerData/config.txt', 'w') as configfile:config.write(configfile)
        configfile.closed
        
##############################################################################################################
    
    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description=""" Starting capture and compression.
        """)
    
    # Add a mutually exclusive for the parser (the arguments in the group can't be given at the same time)
    arg_group = arg_parser.add_mutually_exclusive_group()
    
    arg_group.add_argument('-d', '--duration', metavar='DURATION_HOURS', help="""Start capturing right away, 
        with the given duration in hours. """)
    arg_group.add_argument('-l', '--loops', metavar='NUMBER_LOOPS', help="""Start capturing right away, 
        with the given duration in loops. """)
    
    # Parse the command line arguments
    cml_args = arg_parser.parse_args()
    
    
    # If the duration of capture was given, capture right away for a specified time in hours
    if cml_args.duration:
        
        # Convert the given duration to floating point
        duration = float(cml_args.duration)
    
    elif cml_args.loops:
        
        # Convert the given number of loops to floating point hours by dividing loops by loops per hour
        duration = (float(cml_args.loops))/7.2



    # If a duration was given in loops or hours, run for the desired length of time
    if (cml_args.duration or cml_args.loops):
        
        
        stored_path = makeDirectory(station_code, channel, os.path.join(path, raw) + "/")
        
        # ads1256.run returns 0 when its done recording
        running = ads1256.run(duration, mode, station_code, channel, latitude, longitude, elevation, instrument_string, stored_path)
        
        # Compress all files saved
        source = stored_path.replace("//","/")[:-1]
        
        # Takes the 
        archiveDir(source, os.listdir(source), os.path.join(path, zipped) , os.path.split(source)[1])
        
        if running == 1:
            print("\nProgram terminated!\n")
            sys.exit()
    
    # If no duration given, it is assumed to wait for night and run all night.
    else:
        while(1):
            
            start_time,duration = captureDuration(latitude,longitude,elevation)
            # Convert duration to hours
            duration/=3600
            
            # Wait to start capturing
            if start_time != True:
                
                # Calculate how many seconds to wait until capture starts, and with for that time
                time_now = datetime.datetime.utcnow()
                waiting_time = start_time - time_now
                
                # Wait until sunset
                time.sleep(int(waiting_time.total_seconds()))
    
            path = makeDirectory(station_code, channel, os.path.join(path, raw) + "/")
    
            # ads1256.run returns 0 when its done recording
            running = ads1256.run(duration, mode, station_code, channel, latitude, longitude, elevation, instrument_string, stored_path)
            
            # Compress all files saved
            source = stored_path.replace("//","/")[:-1]
        
            # Takes the 
            archiveDir(source, os.listdir(source), os.path.join(path, zipped) , os.path.split(source)[1])
            
            if running == 1:
                print("\nProgram terminated with ctrl c!\n")
                sys.exit()
            else:
                continue
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
