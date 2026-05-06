# Numerics

This directory contains numerical helpers that are not tied to one physical
subsystem: interpolation, linear solves, convergence tests, excitation helpers,
collision-coefficient interpolation, and statistical-equilibrium assembly.

Keep routines here small and testable. If a routine starts to know about UV
fields, columns, species chemistry, or coolant file layout, it probably belongs
in a physics subdomain instead.

## Background

The model repeatedly solves coupled algebraic and ODE-like problems: molecular
level populations from statistical equilibrium, chemistry abundance evolution,
interpolated collision rates from tabulated molecular data, and convergence
criteria over iterative thermal balance.

Collision interpolation is driven by LAMDA-style molecular data. Statistical
equilibrium and non-LTE line analysis use the same physical inputs as RADEX-like
LVG treatments: Einstein coefficients, collisional coefficients, local density,
radiation field, and escape probability.

## Selected References

- Schoier et al. 2005, LAMDA molecular data, A&A 432, 369,
  https://doi.org/10.1051/0004-6361:20041729
- van der Tak et al. 2007, non-LTE statistical-equilibrium line analysis,
  A&A 468, 627, https://doi.org/10.1051/0004-6361:20066820
- Hindmarsh et al. 2005, SUNDIALS solver suite used by chemistry integration,
  ACM TOMS 31, 363, https://doi.org/10.1145/1089014.1089020
