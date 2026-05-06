# Initialization

This directory prepares the model state before the iterative evolution begins.
It allocates global arrays, constructs spatial indexing, maps particle/grid
data into runtime structures, initializes geometry, and prepares coolant and
chemistry storage.

Initialization code should avoid physical formulas unless they are part of
setting a starting condition. Shared model state lives in `src/physics/state`,
and physics kernels live under `src/physics`.

## Background

3D-PDR uses a point-based representation of the cloud and then casts rays
through the surrounding material to estimate visual extinction, columns, and
radiative coupling. The initialization phase must therefore keep three ideas in
sync:

- point identifiers used by the spatial index,
- PDR/dark-region subsets used by the evolution loop,
- per-ray storage used by column, shielding, and LVG calculations.

## Selected References

- Bisbas et al. 2012, 3D-PDR geometry and point-based PDR modeling, MNRAS 427,
  2100, https://doi.org/10.1111/j.1365-2966.2012.22077.x
- Gorski et al. 2005, HEALPix angular discretization, ApJ 622, 759,
  https://doi.org/10.1086/427976
