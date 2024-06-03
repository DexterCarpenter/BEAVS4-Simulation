# BEAVS 4.0 Simulation 3

Dexter Carpenter

Simulation 3 aims to improve upon simulation 2 by incorporating a basic feedback script into the simulation. The idea is to build a framework of functions that can provide feedback to the system and actively change the blade actuation. A simple function GetFeedback.m is used to calculate a predicted apogee and alter the blade extension based on the difference between the predicted apogee and target apogee. Future iterations will aim to use a PID controller to adjust the blade extension as to avoid undershooting the target apogee.

In addition to a basic feedback script, this simulation also attempts to use Runge Kutta (4,5) via MATLAB's ode45() function. This was done to justify the accuracy of the Fordward Euler method. Results show that the results are extremely similar, which is reassuring. 
