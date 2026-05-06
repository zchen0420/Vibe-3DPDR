# Physics

This directory contains the reusable scientific model pieces. It is split into
state definitions, shared constants, chemistry, radiation, geometry, and
numerical helpers.

The intent is that high-level evolution code can read as a sequence of domain
operations, while the formulas and low-level data layouts stay here.

## Subdirectories

- `core`: physical constants, species indices, global parameters, and the
  central runtime module.
- `state`: derived types for chemistry, coolants, grid, geometry, and thermal
  state.
- `chemistry`: reaction rates, heating rates, H2 formation, and CVODE network
  systems.
- `radiation`: column integration, shielding, dust temperature, photo-rates,
  and LVG escape probabilities.
- `geometry`: HEALPix and ray/evaluation-point helpers.
- `numerics`: interpolation, linear solves, convergence helpers, collision
  interpolation, and statistical-equilibrium helpers.

## Selected References

- Bisbas et al. 2012, 3D-PDR, MNRAS 427, 2100,
  https://doi.org/10.1111/j.1365-2966.2012.22077.x
- Tielens & Hollenbach 1985, PDR chemistry and thermal balance, ApJ 291, 722,
  https://doi.org/10.1086/163111
