# Basic Consepts

## Pin Diode

pin_diode_equivalent_circuit.png<img width="526" alt="image" src="https://user-images.githubusercontent.com/5185118/160285522-6e7d0319-4caf-4839-8c64-f82b7779c21c.png">

Take note of the direction of the current source and the polarity of the diode. Cj depends on reverse voltage of the diode and in this case I have assumed a typical value of 70 pF (data sheet PBW34). For Rs I have asumed 10 ohms, while Rsh is set to 10 Mohm (Mega Ohm).

Both the capacitance of the diode as well as the dark current depend on the reverse bias volatge both of them increasing as the bais volatage increases.

PC mode here means that the cathode is connected to the (inverting) input of the amplifier.
PV mode here means that the anode is connected to the (non_inverting) input of the amplifier.

## Simulation Setup

I have used a poor man's simulator (iCircuit) to demonstrate the basic 2018 circuit using the transimpendance amplifier and the diode in PC mode with 0V bias. In this circuit I have included C1 (4.7 nF) to reduce the bandwidth of the amplifier. ( ~ 1291 Hz based on GBP of 50000 (= GBP of the LMC6464) ).

PC_mode_0V_bias.png<img width="828" alt="image" src="https://user-images.githubusercontent.com/5185118/160283350-a40e22d0-9363-4e13-9305-979f7bc46ca1.png">

The current source is set to 1uA so the expected output is 2.010V. I have included a small ac sweep ontop of the 1uA to study the small signal responce. In this simulation the opamp is considered to be ideal and will not accurately represent the actual system especialy with really small currents and a Vout that is close to zero.

## Simualtion Results

### PC MODE, 0V bias

PC_mode_0V_bias_simulation.png<img width="1440" alt="image" src="https://user-images.githubusercontent.com/5185118/160284559-359a9a8a-d56d-4563-91f6-2942b1765a93.png">

### PV MODE, incorrect !

PV_mode_inverting_simulation.png<img width="1440" alt="image" src="https://user-images.githubusercontent.com/5185118/160285088-1fa944f7-283d-4fc2-98e9-628ee0a4f430.png">

### PV MODE, using dual supply volatge: [+5V, -5V]

PV_mode_dual_supply_simulation.png<img width="1440" alt="image" src="https://user-images.githubusercontent.com/5185118/160284702-d28058d3-e224-4516-b093-bb6379a7f853.png">

### PV MODE, correct ?

PV_mode_non_inverting_simulation.png<img width="1440" alt="image" src="https://user-images.githubusercontent.com/5185118/160284858-19b3e263-0b31-4c44-ac62-015a1b9ea422.png">

To do: 
- Find a better amplifier circuit using a single supply [+5V] that avoids operation near the zero volt supply rail, and the associated non-linearities.

## Test Setup

The current labratory test setup is aimed at comparing various different amplifier and wiring verions. It consist of a aluminiam box housing a raspberry pi 4, a waveshare AD/DA converterboard and a small breadboard, with two identical versions of the 2018 of one of the 2018 amplifier transimpedance amplifiers with 3 diodes each. The ADC is used in differential mode and the V_in_neg is of each differential pair is connected to a Vref that is 1/2 of the summply voltage. The voltage ref IC also provides a temperature output and this used to monitor the temperature of te test box. To prove illumination in a controlled manner there is an led that can be switched on and off. The light of the led is reflected back on the photo diodes via a whitepaper glued on the lid of the box. A simple python script is used to initialize the ADC and read the data, as well as controlling the led.

radiometer_test_setup.png<img width="896" alt="image" src="https://user-images.githubusercontent.com/5185118/160288355-d73ab43f-62c5-45a3-89e2-847c5b3860a0.png">


test_setup_output.png<img width="406" alt="image" src="https://user-images.githubusercontent.com/5185118/160288627-b4c79ee4-c39a-4b9d-86b9-4f7cb065c754.png">

radio_meter_lab_test.png<img width="954" alt="image" src="https://user-images.githubusercontent.com/5185118/160288922-33a9d90f-2454-4f91-a1cf-524a6732278b.png">
