import os
import sys
from collections import OrderedDict

try:
    # Python 2
    import ConfigParser as configparser
except:
    # Python 3
    import configparser




class RDMConfig(object):
    
    def __init__(self):
        self.station_code = None
        self.channel = None
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.instrument_string = None 
        self.raw = None
        self.zipped = None
        self.mode = None
        self.gain = None

        self.upload_enabled = None
        self.hostname = None
        self.rsa_private_key = None
        self.upload_queue_file = None
        self.remote_dir = None
        
        self.read_from_server = None

def readConfig(config_file_path):
    """ Generates two plots of the nights data. 

        Input Arguments:
            
            -config_file_path (string): The path to the directory that stores the configuration file. Ex: /home/pi/RadiometerData/config.txt

        Outputs:

            - rdm_config (object): The configuration object
    """    
    # Create the configuration object
    rdm_config = RDMConfig()
    
    # Create a config object
    config = configparser.ConfigParser()
    
    # Read the config file into the object
    config.read(config_file_path)
    
    # Gather configuration data for the station
    rdm_config.station_code = config['Station']['StationCode']
    rdm_config.channel = config['Station']['Channel']
    rdm_config.latitude = float(config['Station']['Latitude'])
    rdm_config.longitude = float(config['Station']['Longitude'])
    rdm_config.elevation = float(config['Station']['Elevation'])
    rdm_config.instrument_string = config['Station']['InstrumentString']
    rdm_config.raw = config['Station']['RawData']
    rdm_config.zipped = config['Station']['StoredData']
    rdm_config.mode = int(config['Station']['DifferentialMode'])
    rdm_config.gain = int(config['Station']['Gain'])
        
    # Gather configuration data for the upload manager
    rdm_config.upload_enabled = config['Upload']['EnableUpload']
    rdm_config.hostname = config['Upload']['HostName']
    rdm_config.rsa_private_key = config['Upload']['RSAPrivateKey']
    rdm_config.upload_queue_file = config['Upload']['QueueFilename']
    rdm_config.remote_dir = config['Upload']['RemoteDirectory']
    
    rdm_config.read_from_server = config['Server']['ReadFromServer']
    # Return the configuration object
    return rdm_config
        
def makeConfig(config_file_path):
    """ Generates two plots of the nights data. 

        Input Arguments:
            
            -config_file_path (string): The path to the directory that will store the configuration file. Ex: /home/pi/RadiometerData/config.txt

        Outputs:

            - One config.txt file saved in config_file_path
    """    
    # There was no detected config file so one will be created
    # An error message explaining the issue
    print("No config file detected in /home/pi/RadiometerData")
    print("A default config file has been created and can be changed in RadiometerData")
    
    # Create a config object
    config = configparser.ConfigParser()
    
    # optionxform prevents it from naming all config parameters with lower case letters
    config.optionxform = str
    
    # Creates the station data inside the config file using default values
    config['Station'] = OrderedDict((
        ('StationCode', 'AA0000'),
        ('Channel', 'A'),
        ('Latitude', '0.0'),
        ('Longitude', '0.0'),
        ('Elevation', '0.0'),
        ('InstrumentString', 'Your description'),
        ('RawData','CapturedData'),
        ('StoredData','ArchivedData'),
        ('DifferentialMode','1'),
        ('Gain','1')
    ))
        
    # Creates the upload manager configuration section using default settings
    config['Upload'] = OrderedDict((
        ('EnableUpload', 'True'),
        ('HostName', ''),
        ('RSAPrivateKey', '~/.ssh/id_rsa'),
        ('QueueFilename','FILES_TO_UPLOAD.inf'),
        ('RemoteDirectory','.')
    ))

    # Creates the upload manager configuration section using default settings
    config['Server'] = OrderedDict((
        ('ReadFromServer', 'True')
    ))
        
    # Generate the file in the desired directory and close it
    with open(config_file_path, 'w') as configfile:config.write(configfile)
    configfile.closed
    
    # Allow the user to configure the config file
    os.chmod(config_file_path, 0o777)
    
    # Exit allowing the user to configure their settings
    sys.exit()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
