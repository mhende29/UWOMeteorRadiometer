# Notes

These are my notes about the Meteor Radiometer. The first goal was to recreate what was done before: A raspberrypi 4 + WaveAhare high resolution AD and DA converter board + the amplifier based on BPW34 pin diode and LMC6464 quad opamp as described in this paper: 2018_WGN_Radiometer-final.pdf.

WaveShare High-Precision AD/DA Board:

https://www.waveshare.com/high-precision-ad-da-board.htm

https://www.waveshare.com/wiki/High-Precision_AD/DA_Board

The ADC part of the board is based on the ADS1256 from TI:

https://www.ti.com/product/ADS1256

https://www.ti.com/lit/gpn/ads1256

Measuring Single-Ended 0- to 5-V Signals with Differential Delta-Sigma ADCs:

https://www.ti.com/lit/an/sbaa133a/sbaa133a.pdf?ts=1647896597776&ref_url=https%253A%252F%252Fwww.ti.com%252Fproduct%252FADS1256

The amplifier is based on the LMC6464 and the BPW34:

https://www.ti.com/product/LMC6464

https://www.ti.com/lit/gpn/lmc6464

Vishay: https://www.vishay.com/docs/81521/bpw34.pdf

OSRAM: https://dammedia.osram.info/media/resource/hires/osram-dam-5488305/BPW%2034_EN.pdf

General Information on photo diodes:

http://www.osioptoelectronics.com/application-notes/AN-Photodiode-Parameters-Characteristics.pdf

I am now at a stage where I have all of the basic fucntionality working (hardware and software) and I have carried out some basic measurements. Below or some details I have noted down along the way as well as some ideas of what could be done next based on the discussions and suggestions in the [MRM] group as well as my own obersavations.

In order to reduce the noise I have adopted a 2.5V reference source for the Vin_neg of the differential inputs of the ADC instead of a simple resistive voltage divider dirctly from the supply voltage:

https://www.analog.com/media/en/technical-documentation/data-sheets/REF43.pdf

## The Initial Test Setup

RadiometerTestSetup.png<img width="475" alt="image" src="https://user-images.githubusercontent.com/5185118/159028950-c0f64f46-6936-40a8-95d9-bb533a2931c7.png">

## The Amplifier

The amplifiers are setup as a high-sensitivty current to voltage converter.

high_sensitivity_I2V.png<img width="450" alt="I to V converter" src="https://user-images.githubusercontent.com/5185118/159007061-cd312148-de5e-49f8-bb5a-6eac18989d42.png">

    Where: Vout = (-Req) * Iin and Req = (1 + R2/R1 + R2/R3) * R1

    When R1 = 1M, and R2,R3 = 10k this will give a 2.01V/1uA

(Source: "Design with Operational Amplifiers and Analog Intergated Circuits" by Sergio Franco)

The LMC6464 has 4 opamps. In this design 3 of the 4 opamps are setup as high sesnitivty I to V converters each with a 2.01V/1uA output. The 4th opamp is setup as a summing amplifier with a variable gain. The summing amplifier adds the three outputs of the I to V converters and amplifies the signal.

### Possible Improvements

- add a 50Hz | 60z active notch or band reject filter (e.g. Wien-Robinson filter)
- add a anti-aliassing filter (e.g. 4th order Butterworth low pass filter)
- design a new pcb with optimized layout to reduce interference
- adjustable diode bias
- add a 1uA current source to aid with calibration of the amplifier(s)
- instead of using 1/4 of the LMC6464 as a summing amplfier it should also be possible to do the summing of 3 or possibly 4 channels in the software. 
  This should avid the risk of driving the otput of this summing amplifier to the rail voltage and allow for a larger dynamic range of the whole system       withhout sacreficing the single 5V supply setup.

## The Light Sensor

The current source for the intended purpose here is the 3 BPW34 diodes wired in a paralell. The diodes can be used in two ways: PV (Photo Voltaic?) and PC (Photo Conductive?) mode. In this design the diodes are used in the PC mode, meaning that the cathodes will be connceted to the inverting pin of the opamp.
In the group there has been some discussion on whether PC or PV mode is more sensitive. Perhaps a future iteration of the design should have a method for a variable diode bias. Testing with various bias settings could settle this discussion.

The BPW34 DS shows a lower limit of ~0.8uA @ 10lx, while a full moon with a clear sky would between 0.05 to 0.3lx, alterantively a lower limit of ~0.5uA @ 0.01 mW/cm2 and lamba=950nm. (Vishay DS)

### Possible Improvements

- add more diodes (The FDS1010 suugested else where has an active area of 10mmx10mm, while the BPW34 has an active area of 2.62mmx2.65mm. Even with 9 diodes the total area is still 40% less then of the FDS1010.) [I think it is a good idea to stick to individual (smaller) areas as cosmic rays will only pass through one diode and do not affect the whole device.]
- add a tempature sensor and (sofware) compensatation of the temperature dependecie of the BPW34.



## The Analog to Digital Converter

At the moment the ADC is used in a single channel mode. The Vout of the amplifier connects to the Vin_pos of the ADC. This limits the dynamic range of the ADC by 50%. By changing this to the differential mode the full dyanimc range of 24 bits will be enabled, however the number of availble channels will be reduced to 4.

It appears that with the current implemenatation the max smaple rate achievable is ~ 2000 sps, although the ADC is set to a smaple rate of 3750 sps.

### Possible Improvements

- add a Vref source @ 1/2 VDD and bias the Vin_neg @ Vref
- use the DAC in differential mode (software)
- improve the SPI communication with the ADS1256, to improve the maximum achievable sample rate
- add an accurate timing reference (e.g. a GPS 1us pulse) to allow accurate allignment of measurements
