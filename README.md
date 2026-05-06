# Vibe-3DPDR

This repository is a vibe-coding refactor based on my short participation in
the original project. The goal is not to change the scientific model casually,
but to make the code easier to read, build, test, extend, and discuss.

## What This Vibe Edition Adds

The current work focuses on engineering clarity around the existing Fortran
codebase:

- Reorganized the source tree around clear responsibilities: input/output,
  initialization, evolution flow, and physics kernels.
- Split the former flat `physics` directory into domain subdirectories:
  chemistry, core state, geometry, numerics, radiation, and state types.
- Moved global runtime, grid, thermal, geometry, chemistry, species, and coolant
  state into more explicit modules and derived types.
- Separated the main program into a small orchestration entry point at
  `src/main.F90`.
- Added a structured `build/` layout for object files, module files, test
  binaries, and check logs.
- Added focused unit tests for extracted numerical and convergence helpers.
- Tightened many module imports with explicit `only` lists to make dependencies
  easier to understand.
- Renamed legacy abbreviated files and procedures toward lowercase
  `snake_case` names, for example `reaction_rates`, `heating_rates`,
  `collision_coefficients`, and `level_population_system`.
- Replaced the hardest-to-read physics entry points with module APIs that carry
  local state in small derived types: `calculate_collision_coefficients`,
  `calculate_lvg_transition_rates`, `calculate_reaction_rates`, and
  `calculate_heating_rates`.

Readers can find the calculation phase they care about, run tests
before changing behavior, and see where future modularization should happen.

## Source Organization

- `src`
  - `main.F90`: program entry point and top-level orchestration.
  - `io`: parameter handling, input readers, runtime configuration, and final
    output writers.
  - `init`: memory allocation, initial conditions, spatial indexing, particle
    storage, and geometry setup.
  - `evolution`: simulation setup, chemistry refreshes, level population
    solving, thermal balance, and convergence checks.
  - `physics`
    - `core`: shared definitions, global parameters, species indices,
      HEALPix types, and the main state module.
    - `state`: derived types for chemistry, coolant, geometry, thermal,
      and grid state.
    - `chemistry`: reaction rates, heating rates, H2 formation, CVODE
      abundance integration, ODE systems, and Jacobians.
    - `radiation`: columns, radiation field helpers, shielding,
      photo-rate interfaces, dust temperature, and escape probability.
    - `geometry`: HEALPix utilities and evaluation point construction.
    - `numerics`: convergence helpers, excitation helpers, linear
      solves, sorting, splines, and collision coefficient interpolation.
- `tests`: focused unit tests for extracted helpers.

## Build And Use

Requirements:

- `gfortran`
- `gcc`
- SUNDIALS/CVODE installed under `.bin/sundials` by default

Build the executable:

```sh
make
```

Run unit tests:

```sh
make unit-test
```

Run the default regression check:

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

The main compile-time switches live in `Makefile`:

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

## Inheritance And Thanks

This project inherits from the original 3D-PDR code developed by Thomas G.
Bisbas, Tom A. Bell, and collaborators listed in the original program banner.
The scientific foundation belongs to that lineage.

I am especially grateful to Thomas for letting me participate in the project. That experience made this vibe edition possible.

After this Fortran refactor is completed, I plan to create a  Python version
to continue exploring a more usable philosophy of scientific code organization.
Feel free to leave any comment on this project.
