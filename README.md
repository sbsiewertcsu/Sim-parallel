# Sim-parallel
Examples of parallel simulation for planes, trains, automobiles, and spacecraft! The purpose of this code is to provide core integrator/propagator code written for mulit-core and cluster scaling as well as Monte Carlo analysis.

# Train-sim
Starting point for shared memory train-sim - parallel propagator to which Monte Carlo MPI cluster scaling features can be added along with updates to improve fidelity of the simulation and options (e.g., friction, drag, wheel-slip, and other perturbations to the basic acceleration profile for the ideal train that has acceleration from engine thrust - wheel rotation only). The modifications are based upon train simulation theory and practices outlined here - https://www.ecst.csuchico.edu/~sbsiewert/csci551/documents/Lectures/Tutorial-on-Train-Simulation-Parallel-Methods.pdf

Gemini query on MATLAB models for trains (used to verify) - https://www.ecst.csuchico.edu/~sbsiewert/csci551/documents/Papers/train-dynamics/MATLAB-tools-for-simulating-and-analyzing-train-systems.pdf

More references for train simulation - https://www.ecst.csuchico.edu/~sbsiewert/csci551/documents/Papers/train-dynamics/

(CSCI 551 students - note that you can use this code, but it does not solve any exercise questions for MPI simulation of trains - it uses OpenMP for parallel propagation and MPI for Monte Carlo only)
