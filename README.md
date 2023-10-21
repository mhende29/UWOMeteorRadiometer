# UWOMeteorRadiometer

Only supports Python 3, Python 2 not supported. Note: For now, it seems it doesn't work with Pi OS 'Bullseye' nor 'Bookworm'! Please check the wiki for more info.


## Installation

First, you will need to install the BCM2835 library from the Waveshare website, which is used the by ADC: [BCM2835](http://www.airspayce.com/mikem/bcm2835/bcm2835-1.73.tar.gz). 
```
wget https://www.airspayce.com/mikem/bcm2835/bcm2835-1.73.tar.gz
```

Unzip it and follow the instructions in the 'INSTALL' file that should be present in the extracted directory. It's a very basic .configure, make, make install procedure.
```
tar zxvf bcm2835-1.73.tar.gz 
cd bcm2835-1.73
sudo ./configure && sudo make && sudo make check && sudo make install
```

Clone the UWOMeteorRadiometer repository.
```
got clone https://github.com/Habraken/UWOMeteorRadiometer
```

Next, in the UWOMeterRadiometer directory run the ```sudo ./python3_radiometer_installation.sh``` script to install all required packages. This will take a while.

Then in the UWOMeteorRadiometer folder run:
```
sudo python3 setup.py install
```


## Initial setup and recording

After installation, run the radiometer record script:
```
sudo python3 RadiometerRun.py
```
use ```--help```for options.

Sudo privileges are required by the BCM2835 library and there is no way around them.

After the first run the script will tell you that a config file is needed, otherwise upon a second run a default config file will be created. If you don't have a config file, run the script again and it will be created under ~/RadiometerData/config.txt.

Then edit the station code in the config file and the geo coordinates of the radiometer, and the instrument string (i.e. the description of the station).
Then you can run the record script and it will handle the recording automatically.


## Live view

It is also possible to view the signal from the radiometer directly by running:
```
sudo python3 RadiometerLiveView.py
```

This will show an oscilloscope-like window on the screen. You can scroll the mouse inside the window to zoom in the X axis, and press 'A' to automatically scale the Y axis.


## Viewing recorded data

All raw data is stored under ~/RadiometerData/CapturedData as binary RMD files. At the end of every night data is compressed and stored under ArchivedData. Futhermore, 2 plots are generated: a NightPlot file which shows average and maximum levels recorded by the radiometer is a 1 second time interval, and a MaxMinus plot which shows the difference between the maximum level and the average of the neighbouring average levels: max_i - (avg_(i - 1) + avg_(i + 1))/2.
If server upload is configured, data will be uploaded to the server at the end of every night.

The viewing script is not designed to view the RDM files directly, but to query a specific time and the viewing script will extect, join and cut appropriate RDM files and show the signal on the screen.

For example, if there was an event on July 22 2018 at 02:35:21 UTC, and we want to view data from station CA0001 channel A in a 100 second time interval, we would run:
```
python3 AnalyzeData.py CA0001 A 20180722-023521 100
```
