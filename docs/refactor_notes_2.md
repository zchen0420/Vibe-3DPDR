# Refactor Notes 2

This pass scans the larger `src` files after the first split. It focuses on
duplication, naming cleanup, and low-risk consolidation targets rather than
numerical changes.

## Current Pass Update

Completed in this pass:

- Added `src/physics/geometry/ray_path.F90` for shared ray path helpers:
  `ray_step_length`, `projected_point_id`, origin projection initialization,
  projected-point storage, and adaptive-step validation.
- Added `tests/test_ray_path.F90` to lock down the shared indexing and distance
  behavior with a synthetic path.
- Replaced repeated adjacent `epoint` distance formulas in `columns` and
  `radiation` with `ray_step_length`.
- Replaced scattered `projected(:,0)` setup in column, radiation, and level
  population paths with `set_ray_origin_projection`.
- Collapsed the dark-region and PDR-source evaluation point construction into
  one `build_evaluation_points_for_source` helper.
- Avoided duplicate `raytype` allocation in `evaluation_points`; the routine
  now reuses already allocated storage when `particle_storage` owns it.
- Simplified the LVG level-population driver entry points: the coolant solve now
  reads as collision rates, LVG transition rates, then statistical equilibrium,
  instead of passing every cooling-table, ray, and output array at the call site.
- Added `level_population_system_module` with
  `solve_statistical_equilibrium(transition, density, solution)` and removed the
  old `solve_level_population_system` compatibility wrapper.
- Added `point_reaction_rates_module` so chemistry iteration, thermal balance,
  and dark-region chemistry call `calculate_point_reaction_rates(...)` instead
  of each reconstructing the long `calculate_reaction_rates` argument list.
- Added `tests/test_level_population_system.F90` to lock down the new statistical
  equilibrium interface on a small two-level system.
- Replaced the seven coolant collision-partner arrays with
  `collision_rates(partner_id, upper, lower, temperature_index)` plus
  `collision_temperatures(partner_id, temperature_index)`.
- Replaced `READINPUT` with the lowercase module entry
  `read_lamda_coolant_file(coolant_table)`, backed by shared
  `read_collision_partner` and `fill_reverse_collision_rate` helpers.
- Replaced `find_collision_coefficients` with
  `calculate_collision_coefficients(coolant_table, gas_temperature,
  collider_density, coefficients)` and added
  `tests/test_collision_coefficients.F90`.

Updated recommendation:

- Keep generated or network-derived `.c` files out of manual refactor passes.
- Continue with one narrow Fortran domain at a time: shielding naming,
  LAMDA/RADEX collision readers, output writer ordering, chemistry/heating
  branch helpers, then numerics naming.

## Largest Files

The largest files by line count are:

- `src/physics/chemistry/jacobian_full.c` at about 7170 lines
- `src/physics/chemistry/jacobian_mynetwork.c` at about 1367 lines
- `src/physics/chemistry/odes_full.c` at about 701 lines
- `src/physics/chemistry/jacobian_reduced.c` at about 687 lines
- `src/physics/chemistry/calculate_abundances.c` at about 435 lines
- `src/physics/geometry/healpix.F90` at about 398 lines
- `src/physics/chemistry/heating_rates.F90` at about 369 lines
- `src/physics/chemistry/reaction_rates.F90` at about 362 lines
- `src/physics/radiation/h2_shielding.F90` at about 320 lines
- `src/physics/geometry/evaluation_points.F90` at about 294 lines
- `src/physics/numerics/spline.F90` at about 277 lines
- `src/io/output.F90` at about 270 lines
- `src/io/read_input.F90` at about 259 lines
- `src/evolution/level_population_solver.F90` at about 234 lines

The generated or network-derived C files should probably be handled separately
from ordinary refactoring. The best next manual targets are `evaluation_points`,
`h2_shielding`, `co_shielding`, `read_input`, `output`, `heating_rates`, and
`reaction_rates`.

## Priority 1: Shared Ray Path Helpers

`epoint`, `epray`, and `projected` are still manipulated directly in several
places:

- `src/physics/geometry/evaluation_points.F90`
- `src/physics/radiation/columns.F90`
- `src/physics/radiation/radiation.F90`
- `src/evolution/level_population_solver.F90`
- `src/init/geometry_setup.F90`
- `src/init/particle_storage.F90`

Status: mostly complete for path math and indexing. The new `ray_path_module`
now owns shared distance and projected-point helpers:

- `ray_step_length(point, ray_index, eval_index)`
- `projected_point_id(point, ray_index, eval_index)`
- `set_ray_origin_projection(point, point_id)`
- `store_projected_point(...)`
- `update_minimum_adaptive_step(point, adaptive_minimum)`

Remaining work is mostly ownership cleanup: `geometry_setup` and
`particle_storage` still directly allocate or initialize ray storage. Leave that
for a separate memory/state pass so the behavior stays easy to review.

Recommended test coverage:

- Done: add a unit test for `ray_step_length` with a tiny synthetic path.
- Extend the default regression check to assert the printed evaluation point
  count, not only final iteration count.

## Priority 2: Evaluation Point Builder

Status: complete for the duplicated source construction paths.
`src/physics/geometry/evaluation_points.F90` now uses one helper for both:

- the single dark-region source
- the loop over PDR sources

- `build_evaluation_points_for_source(source_point_id)`

That helper now owns:

- allocating and clearing `sorted_distance`, `sorted_point_ids`, and
  per-ray evaluation-point workspace
- filling `sorted_point_ids`
- sorting by distance
- initializing `epoint(:, :, 0)` and `epray`
- walking candidates and stopping rays on ionized points

The top-level routine is now orchestration: optional dark source, all PDR
sources, count/report, negative-step validation, and raytype assignment.

## Priority 3: Shielding Formatting and Names

The first split left `src/physics/radiation/h2_shielding.F90` and
`src/physics/radiation/co_shielding.F90` with inherited indentation from the old
monolithic file. Functions such as `H2SHIELD1`, `H2SHIELD2`, `SCATTER`,
`XLAMBDA`, `COSHIELD`, and `LBAR` are visually nested even though they are
external procedures.

Cleanup plan:

- Normalize all top-level `FUNCTION` and `END` indentation to column 1 or the
  local two-space style.
- Rename public functions to descriptive lowercase names, while keeping
  compatibility wrappers if needed:
  - `H2PDRATE` -> `h2_photodissociation_rate`
  - `H2SHIELD1` -> `federman_h2_self_shielding`
  - `H2SHIELD2` -> `lee_h2_shielding`
  - `SCATTER` -> `dust_scattering_attenuation`
  - `XLAMBDA` -> `extinction_curve_ratio`
  - `COPDRATE` -> `co_photodissociation_rate`
  - `COSHIELD` -> `co_shielding_factor`
  - `LBAR` -> `co_mean_dissociation_wavelength`
  - `CIPDRATE` -> `carbon_photoionization_rate`
  - `SIPDRATE` -> `sulfur_photoionization_rate`
- Consider moving dust scattering helpers out of `h2_shielding.F90` into a
  shared `dust_extinction.F90`, since CO shielding also depends on `SCATTER`.
- Keep table-backed functions close to `shielding_tables.F90`, or wrap table
  interpolation behind helper names so the call sites do not care which table is
  used.

Recommended test coverage:

- Add direct tests for `dust_scattering_attenuation`, `extinction_curve_ratio`,
  and a representative H2/CO shielding value before renaming behavior.

## Priority 4: LAMDA/RADEX Reader

Status: complete for the collision partner layout and reader duplication.
`src/io/read_input.F90` now reads LAMDA/RADEX collision partners with one helper
instead of seven partner-specific branches. The old repeated sequence was:

- skip a line
- read temperatures
- read collision rows
- store coefficients into one partner-specific array
- calculate reverse rates by detailed balance

This is now handled by helpers:

- `read_collision_partner(input_unit, partner_id, collision_count, temperature_count, coolant_table)`
- `fill_reverse_collision_rate(coolant_table, partner_id, upper_level, lower_level, temperature_index)`

The separate arrays (`H_COL`, `HP_COL`, `EL_COL`, `HE_COL`, `H2_COL`,
`PH2_COL`, `OH2_COL`) were removed from `coolant_data` and replaced by:

- `collision_temperatures(partner_id, temperature_index)`
- `collision_rates(partner_id, upper, lower, temperature_index)`

`src/physics/numerics/collision_coefficients.F90` is now a module exposing
`calculate_collision_coefficients`. It receives a coolant table and local
collider-density vector, so the LVG solver no longer knows about file-layout
details or seven partner-specific arrays.

## Priority 5: Output Writers

`src/io/output.F90` still has repeated file-open and row-writing patterns for:

- `.pdr.fin`
- `.ion.fin`
- `.mol.fin`
- `.cool.fin`
- `.heat.fin`
- `.line.fin`
- `.spop.fin`
- `.opdp.fin`

Good next helpers:

- `open_output_file(unit, suffix)`
- `write_region_rows(unit, region_ids, region_count, include_dust_temperature)`
- `write_pdr_rows_in_boundary_order(unit, row_writer)`
- `write_coolant_series_row(...)`

The special PDR ordering `1:pdr_ptot-2`, then `pdr_ptot`, then `pdr_ptot-1`
appears more than once and should be named. A helper like
`ordered_pdr_point_index(sequence_index)` would make the boundary handling
explicit and testable.

## Priority 6: Heating and Reaction Kernels

The first pass extracted a few kernels, but both main routines remain long.

For `src/physics/chemistry/heating_rates.F90`, the next splits should be:

- dust photoelectric heating into a helper that owns the Newton solve
- PAH heating/cooling into one helper returning net PAH heating plus optional
  diagnostics
- H2 formation, photodissociation, and FUV pumping heating into helpers
- chemical heating per network into `chemical_heating_reduced`,
  `chemical_heating_full`, and `chemical_heating_mynetwork`

For `src/physics/chemistry/reaction_rates.F90`, the next splits should be:

- `apply_photoreaction_rates`
- `apply_duplicate_photorate`
- `apply_cosmic_ray_ionization_rate`
- `apply_cosmic_ray_photon_rate`
- `apply_surface_and_desorption_rates`

The current `GOTO` labels are a good migration target. Replace them with a
small reaction-class dispatch once each branch has a helper.

## Priority 7: Numerics Cleanup

`src/physics/numerics/spline.F90`, `linear_solver.F90`, and
`level_population_diagnostics.F90` still use old uppercase procedure names and
comment styles. Recommended cleanup:

- Partial progress: `level_population_system.F90` now exposes only the
  descriptive module entry `solve_statistical_equilibrium(...)`; the old
  external name was removed.
- Rename `GAUSS_JORDAN` to `solve_linear_system_gauss_jordan`.
- Rename `GAUSS_JORDAN_writes` to `dump_level_population_solver_diagnostics`.
- Rename `SPLINE`, `SPLINT`, `SPLIE2`, and `SPLIN2` to lowercase descriptive
  names, with wrappers if external call sites still use legacy names.
- Replace `call_writes` with a clearer status result such as `singular_matrix`.
- Add a small linear solver unit test before renaming.

## Priority 8: Generated Chemistry C Files

The largest files are chemistry network C files, especially
`jacobian_full.c`. Before hand-editing these, confirm whether they are generated
from reaction-network definitions. If they are generated, prefer:

- documenting the generator and source network files
- excluding generated files from style refactors
- adding checks that generated output is current

If they are not generated, they need a separate plan because the blast radius is
larger than ordinary Fortran cleanup.

## Suggested Order

1. Done: add ray path helpers and tests around `epoint` distance/index behavior.
2. Done: collapse `evaluation_points` dark/PDR construction into one source builder.
3. Extend the default regression check to assert the printed evaluation point count.
4. Normalize shielding indentation, then add wrappers for descriptive names.
5. Done: extract LAMDA collision-partner reader helpers and move the seven
   collision-partner arrays into one `collision_rates(...)` table.
6. Table-drive output writers and name the special PDR output order.
7. Continue shrinking `heating_rates` and `reaction_rates` branch by branch.
8. Clean numerics names after adding focused tests.
9. Keep generated chemistry C files in a separate generator/documentation track.
