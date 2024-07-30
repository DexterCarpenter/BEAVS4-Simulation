# BEAVS 4.0 Simulation

Dexter Carpenter

Winter & Spring, 2024

This repository hosts MATLAB scripts and data to run a simulation for Oregon State University's Experimental Sounding Rocket Association Team. The simulation is used to gain understanding of how the Blade Extending Apogee Variance System (BEAVS) will affect the rocket's flight path. The simulation takes [OpenRocket](https://openrocket.info/) data export in conjunction with [OpenMotor](https://github.com/reilleya/openMotor) custom motors as a basis.

Objectives
- Understand the total braking power of BEAVS given varying rocket designs
- Separation failure risk analysis
- PID control simulation
- Output data to BEAVS hardware to simulate flight (did not achieve)

### [Simulation 1](https://github.com/DexterCarpenter/BEAVS4-Simulation/tree/main/%5B5.8.5.1%5DBEAVS4_Sim1_R2020b#beavs-40-simulation-1)

Preliminary simulation and proof of concept.

### [Simulation 2](https://github.com/DexterCarpenter/BEAVS4-Simulation/tree/main/%5B5.8.5.2%5DBEAVS4_Sim2_R2020b#beavs-40-simulation-2)

Simulation 2 aims to provide a better representation for the Coefficient of Drag. Future iterations of the simulation will aim to control the BEAVS extension.

### [Simulation 3](https://github.com/DexterCarpenter/BEAVS4-Simulation/tree/main/%5B5.8.5.3%5DBEAVS4_Sim3_R2020b#beavs-40-simulation-3)

Incorporates a basic feedback script into the simulation. A function is used to calculate a predicted apogee and alter the blade extension based on the difference between the predicted apogee and target apogee. This version also attempts to use Runge-Kutta (4,5) to solve the ODE.

### [Simulation 4](https://github.com/DexterCarpenter/BEAVS4-Simulation/tree/main/%5B5.8.5.4%5DBEAVS4_Sim4_R2020b#beavs-40-simulation-4)

Adds PID option to start finding Kp, Ki, and Kd values. This PID evaluated the error between the current _altitude_ and a target _altitude_ of 10,000ft.

### [Simulation 5](https://github.com/DexterCarpenter/BEAVS4-Simulation/tree/main/%5B5.8.5.5%5DBEAVS4_Sim5_R2020b#beavs-40-simulation-5)

Instead of having the PID evaluate the error between the _altitude_, it now uses a lookup table to determine the error between a _velocity_. The lookup table is based on a previous test flight.

### Use

To use this simulation to test your own [OpenRocket](https://openrocket.info/) rockets, you must use the following file export format:
> OpenRocket Data Export config:
> * Must go Edit > Preferences > Units > Default Metric > Ok
> * Select All Variables
> * Format > Field Seperator String > Comma
> * Decimal Places > 3
> * Use Exponential Notation
> * Section Comments > Don't Include Simulation Description
> * Section Comments > Don't Include Flight Events
> * Comment Character > #
> * Keep Field Descriptions
    
> OpenRocket Event Export config:
> * Variables to Export > None
> * Format Settings > csv
> * Include Flight Events
> * Comment Character > #

However, this script was not written for general use. The primary function is to provide the OSU ESRA team insights. Modification of the data display requires edits to the script, and is not made easy to use for an outside user. Use outside of the OSU ESRA team is certainly encouraged.

### Reference

See also Colin Hale-Brown's repository for [BEAVS 4.0](https://github.com/colinhalebrown/BEAVS4)

This uses [MATLAB R2020b](https://www.mathworks.com/downloads/)
