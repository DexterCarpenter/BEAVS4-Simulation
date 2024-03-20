# BEAVS 4.0 Simulation

Dexter Carpenter

Winter & Spring, 2024

This repository hosts MATLAB scripts and data to run a simulation for the Oregon State University's Experimental Sounding Rocket Association Team. The simulation is used to gain understanding of how the Blade Extending Apogee Variance System (BEAVS) will affect the rocket's flight path. The simulation takes [OpenRocket](https://openrocket.info/) data export in conjunction with [OpenMotor](https://github.com/reilleya/openMotor) custom motors as a basis.

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

This uses [MATLAB R2020b](https://www.mathworks.com/downloads/)
