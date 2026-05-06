# Radiation

This directory contains radiation-field propagation, column-density integration,
photo-rate interfaces, H2 and CO shielding, dust temperature estimates, atomic
photo-rates, and LVG escape-probability calculations.

## Background

In a PDR, far-ultraviolet photons regulate molecular dissociation, carbon and
sulfur ionization, grain/PAH heating, and therefore the cooling-line emission
that observers use as diagnostics. The code tracks directional columns and
visual extinction along HEALPix rays, applies shielding functions, and uses
escape probabilities to connect local level populations with emergent line
cooling.

The LVG path here is deliberately local: it receives a coolant table, local
conditions, ray samples, level populations, and collision coefficients, then
updates radiative transition rates, optical depths, escape probabilities, line
emission, and cooling.

## Selected References

- Draine & Bertoldi 1996, H2 shielding in photodissociation fronts, ApJ 468,
  269, https://doi.org/10.1086/177689
- Visser, van Dishoeck & Black 2009, CO photodissociation and shielding,
  A&A 503, 323, https://doi.org/10.1051/0004-6361/200912129
- Tielens & Hollenbach 1985, PDR radiation, chemistry, and line diagnostics,
  ApJ 291, 722, https://doi.org/10.1086/163111
- van der Tak et al. 2007, LVG/non-LTE line analysis context, A&A 468, 627,
  https://doi.org/10.1051/0004-6361:20066820
- Gorski et al. 2005, HEALPix directional sampling, ApJ 622, 759,
  https://doi.org/10.1086/427976
