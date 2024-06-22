# BEAVS 4.0 Simulation 5

Dexter Carpenter

Simulation 5 improves upon simulation 4's PID. Instead of the PID error function evaluating the error between the rocket's current altitude and a target altitude of 10,000ft, the simulation 5 error function evaluates the error between the current velocity of the rocket and a desired velocity curve. This velocity curve acts as a lookup table that maps an altitude to a desired velocity, starting at around where motor cutoff should ocurr and ending at an altitude slightly above 10,000ft with a velocity of zero. The PID function in this simulation evaluates this error. The PID constants are tuned for this system.

Furthermore, the figures have been altered to better show how PID is acting over time.
