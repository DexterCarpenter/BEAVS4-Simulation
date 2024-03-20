# BEAVS 4.0 Simulation 2

Dexter Carpenter

This simulation expands upon simulation 1. Simulation 2 aims to provide a better representation for the Cd (coefficient of Drag) of the rocket and BEAVS. Further research has shown that in order to estimate a Cd for the entire rocket, the following equation is used to add the coefficient of drag from BEAVS to the coefficient of drag of the rocket.

`Cd_total = Cd_rocket + Cd_BEAVS * ( A_BEAVS / A_Ref )`

Where `Cd_rocket` is the total Cd from OpenRocket, `Cd_BEAVS` is estimated to be ~1.10 as a flat plate against a wall, `A_BEAVS` is the area of all BEAVS blades combined excluding the body tube, and the `A_ref` is the reference area of the rocket defined by the axiala cross-sectional area at the base of the nosecone.

Like simulation 1, this simulation assumes a maximum extension of BEAVS. However, the extension begins at 0 at the beginning of the coast phase and extends at a masimum rate until fully extended.

Future iterations of the simulation will aim to control the BEAVS extension to achieve a target apogee of 10,000 ft.
