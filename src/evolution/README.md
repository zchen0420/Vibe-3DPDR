# Evolution

This directory owns the iteration-level workflow. It coordinates chemistry
updates, ray-column refreshes, LVG level-population solves, heating/cooling
balance, dark-region handling, and convergence tests.

The files here should read like orchestration. Physics formulas belong in
`src/physics`, while this layer decides when to call them and how to move data
between the grid state, chemistry state, thermal state, and coolant iteration
buffers.

## Background

PDR models are coupled nonlinear systems: chemistry depends on radiation and
temperature, cooling depends on level populations, and thermal balance feeds
back into chemistry. 3D-PDR follows this by iterating chemistry, radiation
columns, statistical equilibrium, cooling, and heating until the model reaches
a convergence criterion.

The LVG work path combines local collision rates with line escape probabilities
and then solves statistical equilibrium for coolant populations. The thermal
path compares total heating against total cooling and updates the gas
temperature through guarded iteration and bisection-like moves.

## Selected References

- Bisbas et al. 2012, 3D-PDR, MNRAS 427, 2100,
  https://doi.org/10.1111/j.1365-2966.2012.22077.x
- van der Tak et al. 2007, "A computer program for fast non-LTE analysis of
  interstellar line spectra", A&A 468, 627,
  https://doi.org/10.1051/0004-6361:20066820
- Tielens & Hollenbach 1985, PDR thermal and chemical balance, ApJ 291, 722,
  https://doi.org/10.1086/163111
