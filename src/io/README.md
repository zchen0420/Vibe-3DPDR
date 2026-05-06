# Input And Output

This directory owns data ingestion and final product writing. It reads runtime
parameters, species and reaction files, LAMDA/RADEX-style coolant data, initial
grid inputs, and writes PDR, cooling, heating, line, level-population, and
optical-depth outputs.

Parsing code should keep file-format knowledge here. Once data is read, callers
should work with typed state such as `reaction_network`, `coolant_data`, and
`runtime_config_state` rather than raw input columns.

## Background

The chemistry network follows the common astrochemical pattern of temperature
bounded Arrhenius-like rates, photoreaction rates, cosmic-ray reactions, and
grain-related reactions. Coolant files use molecular spectroscopy data:
energies, statistical weights, Einstein coefficients, frequencies, and
collisional rate coefficients.

Output files preserve the traditional 3D-PDR product split so results remain
recognizable to existing workflows.

## Selected References

- McElroy et al. 2013, "The UMIST database for astrochemistry 2012", A&A 550,
  A36, https://doi.org/10.1051/0004-6361/201220465
- Schoier et al. 2005, "An atomic and molecular database for analysis of
  submillimetre line observations", A&A 432, 369,
  https://doi.org/10.1051/0004-6361:20041729
- van der Tak et al. 2007, RADEX non-LTE line analysis, A&A 468, 627,
  https://doi.org/10.1051/0004-6361:20066820
