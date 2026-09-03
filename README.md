# BiomolecularCondensateNucleationMC

A Monte Carlo simulation framework for studying the nucleation of biomolecular condensates in the grand canonical ensemble.

## Features

The current code includes:
- Grand canonical Monte Carlo (GCMC)
- Aggregation-volume-bias Monte Carlo (AVBMC)
- Configurational-bias Monte Carlo (CBMC)
- Translation and rotation Monte Carlo moves
- Long-chain insertion and deletion moves
- Umbrella sampling
- Weighted histogram analysis method (WHAM)
- Distance- and energy-based cluster criteria
- Multiple interaction models and force fields

## Interaction Models

Implemented interaction models include:
- Mpipi
- HPS variants
- Lennard-Jones and electrostatic interactions
- Pedone potential
- Tersoff potential

## Source Code

The main Fortran source code is located in the `src/` directory.

## Compilation

The code is written primarily in Fortran and uses MPI compiler wrappers.

To compile the program:

```bash
make
```

The default Makefile produces the executable `generalNucleation`.

To clean generated build files:

```bash
make clean
```

## Input and Simulation Setup

Simulation behavior is controlled through input scripts. Information about available script commands and analysis variables is provided in the `docs/` directory.

## Code Provenance

This project is derived from the `NucleationSimulationMC` code originally developed by Troy D. Loeffler:
https://github.com/mrnucleation/NucleationSimulationMC

The present repository contains subsequent development and modifications for simulations of biomolecular condensate nucleation, including additional Monte Carlo methods, interaction models, and analysis capabilities.

## Authors and Contributors

Current development is carried out in the Zerze research group at the University of Houston.

The original `NucleationSimulationMC` code was developed by Troy D. Loeffler. See the original repository linked above for the code from which this project was derived.
