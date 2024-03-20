# BEAVS 4.0 Simulation 1

Dexter Carpenter

This simulation is a preliminary simulation. It serves as a proof of concept that importing data from OpenRocket is feasible, and data is useable. This simulation implements Forward Euler to calculate a flight path using initial conditions from the OpenRocket data. This simulation assumes the BEAVS system is aucuating at a maximum instantaously from the beginning of the coast phase. This provides a very rough estimate of the theoretical maximum BEAVS could lower the rocket's apogee.

There are many concerns regarding the accuracy of BEAVS effect in this simulation. The Coefficient of Drag is extremely hard to determine without empirical data. This is a primary problem to solve in all simulations. 
