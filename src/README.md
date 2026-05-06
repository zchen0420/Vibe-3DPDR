# Source Tree

This directory contains the Fortran and generated C implementation of the
3D-PDR model. The code is organized around the runtime flow rather than around
the old flat source layout:

- `io`: reads runtime parameters, chemical networks, coolant data, and writes
  model products.
- `init`: allocates state, builds geometry, and prepares point storage.
- `evolution`: advances chemistry, level populations, cooling, heating, and
  convergence decisions.
- `physics`: contains the reusable model kernels, state types, constants, and
  numerical helpers.

The code models photodissociation regions (PDRs): neutral interstellar gas
where far-ultraviolet radiation controls chemistry, heating, cooling, and line
emission. The implementation inherits its scientific framing from the original
3D-PDR code and keeps the refactor intentionally close to that model.

## Style Notes

Fortran source in this tree uses lowercase code, two-space indentation, and
lowercase `snake_case` file names. The generated chemistry C files are kept
separate from manual style work because they represent network-generated
systems rather than hand-maintained application logic.

## Selected References

- Bisbas et al. 2012, "3D-PDR: a new three-dimensional astrochemistry code for
  treating photodissociation regions", MNRAS 427, 2100,
  https://doi.org/10.1111/j.1365-2966.2012.22077.x
- Tielens & Hollenbach 1985, "Photodissociation regions. I. Basic model",
  ApJ 291, 722, https://doi.org/10.1086/163111
- Gorski et al. 2005, "HEALPix: A Framework for High-Resolution Discretization
  and Fast Analysis of Data Distributed on the Sphere", ApJ 622, 759,
  https://doi.org/10.1086/427976
