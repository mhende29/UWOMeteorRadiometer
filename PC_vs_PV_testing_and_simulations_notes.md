# Basic Consepts

## Pin Diode

pin_diode_equivalent_circuit.png<img width="468" alt="image" src="https://user-images.githubusercontent.com/5185118/160282593-c18123a4-d747-49b8-9919-2f52e87aa308.png">

Take note of the direction of the current source and the polarity of the diode. Cj depends on reverse voltage of the diode and in this case I have assumed a typical value of 70 pF (data sheet PBW34). For Rs I have asumed 10 ohms, while Rsh is set to 10 Mohm (Mega Ohm).

Both the capacitance of the diode as well as the dark current depend on the reverse bias volatge both of them increasing as the bais volatage increases.

## Simulation Setup

I have used a poor man's simulator to demonstrate the basic 2018 circuit using the transimpendance amplifier and the diode in PC mode with 0V bias. In this circuit I have included C1 (4.7 nF) to reduce the bandwidth of the amplifier. ( ~ 1291 Hz based on GBP of 50000 (= GBP of the LMC6464) ).

PC_mode_0V_bias.png<img width="828" alt="image" src="https://user-images.githubusercontent.com/5185118/160283350-a40e22d0-9363-4e13-9305-979f7bc46ca1.png">

The current source is set to 1uA so the expected output is 2.010V. I have included a small ac sweep ontop of the 1uA to study the small signal responce. In this simulation the opamp is considered to be ideal and will not accurately represent the actual system especialy with really small currents and a Vout that is close to zero.
