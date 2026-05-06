# State Types

This directory defines the derived types that hold model state: chemistry
networks, coolant tables and per-point coolant state, geometry columns,
simulation-grid points, and thermal arrays.

State types should describe ownership and shape. They should avoid hidden
physics behavior; calculations belong in the domain modules that consume these
types.

## Background

PDR calculations share the same physical state across several coupled solvers:
chemistry needs abundance and local radiation, heating needs reactions and gas
conditions, LVG cooling needs coolant data and level populations, and output
needs all of those fields in stable shapes. Keeping these shapes explicit makes
the coupling visible and keeps long argument lists out of the physics kernels.

## Selected References

- Bisbas et al. 2012, 3D-PDR model state and coupled workflow, MNRAS 427,
  2100, https://doi.org/10.1111/j.1365-2966.2012.22077.x
- McElroy et al. 2013, UMIST reaction-network context, A&A 550, A36,
  https://doi.org/10.1051/0004-6361/201220465
- Schoier et al. 2005, LAMDA coolant data context, A&A 432, 369,
  https://doi.org/10.1051/0004-6361:20041729
