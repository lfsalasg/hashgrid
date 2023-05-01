# Introduction

This example provides a very simple implementation of the Monte Carlo algorithm for the calculation of the energy of a Lennard-Jones fluid in a NVT ensemble. The example uses the 6-12 lennard-jones to represent the van der Waals interactions and the metropolis criteria to accept/reject the movement of a particle in the system.

For a deeper explanation of the monte carlo technique for thermodynamics see [ref](https://en.wikipedia.org/wiki/Monte_Carlo_method_in_statistical_mechanics)

# System representation

The system consists of a cubic box of 27 A with a cutoff of 9 A.  The system is initialized with 100 particles at random position. Then each particle is deffined as a lennard-jones sphere that interacts with each other using a classical 12-6 potential (the parameters are fixed at 158.5 K and 3.72 A)

**Note**: The temperature, box size, cutoff, number of particles and steps can be changed. Run `./target/release/examples/1_montecarlo --help` to see the CLI options.