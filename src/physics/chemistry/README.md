# Chemistry

This directory contains the chemical and heating-rate kernels. It computes
reaction coefficients from a loaded network, evaluates local heating channels,
handles H2 formation rates, and contains generated C systems for CVODE
abundance integration.

The hand-maintained Fortran files are formatted with the project style. The
large `odes_*.c` and `jacobian_*.c` files are generated or network-derived
chemistry systems; treat them as generated scientific artifacts unless the
network generation path is being changed intentionally.

## Background

The chemistry model combines gas-phase reactions, photoreactions, cosmic-ray
processes, freeze-out/desorption terms, grain-surface H2 formation, and
network-specific ODE/Jacobian systems. Heating terms include photoelectric
heating, PAH heating/cooling, carbon ionization, H2 formation and
photodissociation heating, FUV pumping, cosmic rays, turbulent heating,
chemical heating, and gas-grain exchange.

## Selected References

- McElroy et al. 2013, UMIST Database for Astrochemistry 2012, A&A 550, A36,
  https://doi.org/10.1051/0004-6361/201220465
- Tielens & Hollenbach 1985, foundational dense PDR chemistry and heating,
  ApJ 291, 722, https://doi.org/10.1086/163111
- Bakes & Tielens 1994, photoelectric heating by small grains and PAHs,
  ApJ 427, 822, https://doi.org/10.1086/174188
- Weingartner & Draine 2001, grain charging and gas heating, ApJS 134, 263,
  https://doi.org/10.1086/320852
- Wolfire et al. 2003, neutral atomic phases and heating/cooling balance,
  ApJ 587, 278, https://doi.org/10.1086/368016
- Hindmarsh et al. 2005, SUNDIALS nonlinear and differential/algebraic solvers,
  ACM TOMS 31, 363, https://doi.org/10.1145/1089014.1089020
