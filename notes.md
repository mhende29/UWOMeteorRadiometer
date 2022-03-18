# NOTES
These are my notes about the Meteor Radiometer. The first goal was to recreate what was done before: A raspberrypi 4 + WaveAhare high resolution AD and DA converter board + the amplifier based on BPW34 pin diode and LMC6464 quad opamp as described in this paper: 2018_WGN_Radiometer-final.pdf.

WaveShare High-Precision AD/DA Board:

https://www.waveshare.com/high-precision-ad-da-board.htm

https://www.waveshare.com/wiki/High-Precision_AD/DA_Board

The ADC part of the board is based on the ADS1256 from TI:

https://www.ti.com/product/ADS1256

https://www.ti.com/lit/gpn/ads1256

The amplifier is based on the LMC6464 and the BPW34:

https://www.ti.com/product/LMC6464

https://www.ti.com/lit/gpn/lmc6464

Vishay: https://www.vishay.com/docs/81521/bpw34.pdf

OSRAM: https://dammedia.osram.info/media/resource/hires/osram-dam-5488305/BPW%2034_EN.pdf

The amplifiers are setup as a high-sensitivty current to voltage converter.

high_sensitivity_I2V.png<img width="450" alt="I to V converter" src="https://user-images.githubusercontent.com/5185118/159007061-cd312148-de5e-49f8-bb5a-6eac18989d42.png">

Where: Vout = (-Req) * Iin and Req = (1 + R2/R1 + R2/R3) * R1

When R1 = 1M, and R2,R3 = 10k this will give a 2.01V/1uA

(Source: "Design with Operational Amplifiers and Analog Intergated Circuits" by Sergio Franco)

The current source for the intended purpose here is the 3 BPW34 diodes wired in a paralell. The diodes can be used in two ways: PV (Photo Voltaic?) and PC (Photo Conductive?) mode. In this design the diodes are used in the PC mode, meaning that the cathodes will be connceted to the inverting pin of the opamp.
In the group there has been some discussion on whether PC or PV mode is more sensitive. Perhaps a future iteration of the design should have a method for a variable diode bias. Testing with various bias settings could settle this discussion.

The BPW34 DS shows a lower limit of ~0.8uA @ 10lx, while a full moon with a clear sky would between 0.05 to 0.3 lx.
