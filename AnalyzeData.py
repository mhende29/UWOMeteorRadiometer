

from __future__ import print_function

import argparse
import numpy as np


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


def readBinary(file_name, checksum_check=False, read_data=True):
    """ Read a binary file. """

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

        




if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    
    # Set up input arguments
    arg_p = argparse.ArgumentParser(description="Analyzes radiometer files.")
    arg_p.add_argument('rdm_file', type=str, help="Path to the .rdm file.")
    
    # Parse input arguments
    cml_args = arg_p.parse_args()
    
    
    # Read the binary RDM file
    rdm , chksum_pass = readBinary(cml_args.rdm_file)
    
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
    plt.plot(unix_times, rdm.intensity)
    #plt.xlim((1526482266,1526482272))
    #plt.ylim((2740000,2760000))
    print(rdm.num_samples/(unix_times[-1]-unix_times[0]))
    plt.show()
    
    

    
    
    
    
