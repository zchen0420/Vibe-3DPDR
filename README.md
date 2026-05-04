# 3D-PDR

3D-PDR is a photodissociation region simulation code for evolving chemistry,
line level populations, radiative quantities, and thermal balance on a grid.
This repository is now reorganized around three top-level program phases:
reading inputs, initialising memory and calculation structures, and evolving
the simulation to convergence.

This project is fully built via vibe coding based on my short participation in the project.

## Build

Requirements:

- `gfortran`
- `gcc`
- SUNDIALS/CVODE installed under `.bin/sundials` by default

Build the executable:

```sh
make
```

Run the verification suite:

```sh
make check
```

All intermediate files are written under `build/`:

- `build/obj`: object files
- `build/mod`: Fortran module files
- `build/tests`: unit test executables
- `build/logs`: check logs

Clean generated files:

```sh
make clean
```

The `.bin/` directory is intentionally left alone by cleanup because it may
contain the local SUNDIALS runtime.

## Configuration

The default runtime configuration is:

```text
configs/default.params
```

Run with the default configuration:

```sh
./3DPDR
```

Run with another parameter file:

```sh
./3DPDR path/to/your.params
```

The main compile-time switches live in `makefile`:

- `DIMENSIONS`: `1`, `2`, or full 3D mode by omission of pseudo flags
- `NETWORK`: `REDUCED`, `FULL`, or `MYNETWORK`
- `DUST`: dust treatment mode
- `GUESS_TEMP`, `THERMALBALANCE`, `TEMP_FIX`, `CO_FIX`, `H2FORM`: feature flags
- `SUNDIALS_PREFIX`: location of the SUNDIALS installation, default `.bin/sundials`

Example:

```sh
make clean
make NETWORK=FULL DIMENSIONS=1 check
```

## Usage

The default example uses `data/1Dn30.dat` and writes outputs with the `V1`
prefix:

- `V1.pdr.fin`: PDR abundances and temperatures
- `V1.cool.fin`: cooling rates
- `V1.heat.fin`: heating rates
- `V1.line.fin`: transition line outputs
- `V1.spop.fin`: level populations
- `V1.opdp.fin`: optical depths in pseudo-1D mode

The regression check asserts that the default run converges with:

```text
RESULT status=converged iterations=227
```

## Source Layout

- `src/app/3DPDR.F90`: main program and top-level orchestration of the read, init, and evolve phases
- `src/io`: parameter handling, input file readers, and final output writers
- `src/init`: memory allocation, initial conditions, spatial indexing, and geometry setup
- `src/evolution`: iteration setup, chemistry refreshes, level populations, thermal balance, and convergence
- `src/physics`: physical constants, global state modules, radiative/column/cooling helpers, HEALPix utilities, and numerical kernels
- `src/physics/chemistry/*.c`: CVODE abundance integration, ODE systems, and Jacobians
- `tests/*.F90`: focused unit tests for extracted numerical helpers

## Acknowledgements

Original 3D-PDR development is credited to Thomas G. Bisbas, Tom A. Bell, and
collaborators listed in the original program banner. Special thanks to Thomas
for making this codebase possible and helping me.

This refactor was completed by Zhousi Chen through vibe coding, finishing
engineering cleanup and modularisation work that had not been completed before.

Future plans include migrating the workflow toward Python-facing interfaces and
GPU computation.
