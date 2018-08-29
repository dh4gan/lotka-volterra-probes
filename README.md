# Predator-Prey populations of self-replicating interstellar probes evolving on a stellar network

This code (written in C++ using the Eclipse CDT) follows the growth of self-replicating interstellar probes on a network of stars.

Each star solves the Lotka-Volterra equations, where each system is coupled to neighbouring systems by inflow/outflow terms for both prey and predators.  More information can be found in the accompanying paper:

<link to follow>

## Compiling the code ##

This code was written in the Eclipse CDT, and hence requires this to produce Makefiles for compilation.  A standard g++ Makefile is coming!

## Running the code

The code is run using the command

`> ./lotka_volterra_probes <paramfile>`

An example can be found in the `paramfiles/` directory.  The parameters that define each Lotka-Volterra system can be fixed as constants, sampled from a uniform distribution or sampled from a Gaussian.  The parameter file explains how each sampling method is defined/input into the code.

