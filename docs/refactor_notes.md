# Refactor Notes

## Naming Conventions

- Source file names use lowercase `snake_case`.
- Fortran module names use lowercase `snake_case_module`.
- Derived types use lowercase `snake_case`.
- New procedures should prefer descriptive verb phrases, for example `calculate_reaction_rates` rather than abbreviated names.
- Domain abbreviations are kept when they are standard in the model, for example `h2`, `pdr`, `uv`, `lvg`, and `co`.

## Completed Split Work

### Radiation Shielding

The old `src/physics/radiation/shield.F90` routines were split by domain:

- `h2_shielding.F90`
- `co_shielding.F90`
- `atomic_photo_rates.F90`

`shielding_tables.F90` remains shared read-only tabulated data.

### Reaction Rates

`src/physics/chemistry/reaction_rates.F90` keeps the public `calculate_reaction_rates`
entry point, while reusable rate kernels now live in `reaction_rate_kernels_module`.
The extracted helpers cover:

- Arrhenius rate selection
- freeze-out rates
- cosmic-ray and desorption rates
- rate clamping and warning output

### Heating Rates

`src/physics/chemistry/heating_rates.F90` now uses `heating_rate_kernels_module`
for independent terms and the `heating_rate(1:12)` output layout. The output
contract is covered by `tests/test_heating_rate_layout.F90`.

Extracted helpers cover:

- cosmic ray and turbulent heating
- gas-grain exchange
- output layout storage

### Level Population Solver

`src/evolution/level_population_solver.F90` now has a `coolant_work_item` type
that carries:

- coolant id
- species abundance index
- level count
- temperature count
- workspace arrays

The repeated CII, CI, OI, and CO blocks now share `solve_coolant_population`.

### Evaluation Points

`src/physics/geometry/evaluation_points.F90` now shares the ray projection and
projected-point storage path between dark-region and PDR-region construction.
The full default regression check covers projected point behavior.

### Level Population System

`src/physics/numerics/level_population_system.F90` was split into:

- public level population solve wrapper
- `linear_solver.F90`
- `level_population_diagnostics.F90`

The solver should eventually be replaced or isolated behind a small interface, but that is a numerical decision rather than a naming cleanup.
