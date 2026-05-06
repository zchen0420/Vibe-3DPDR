# Core

This directory defines shared constants, global parameters, species index maps,
HEALPix-related scalar types, and the central `maincode_module` state handles.

Core files should stay boring and explicit. They are allowed to know about the
model-wide state layout, but they should not own chemistry, radiation, or
thermal formulas.

## Background

PDR calculations move repeatedly between physical units, species indices,
coolant indices, radiation geometry, and global runtime switches. This layer
keeps those conventions centralized so physics kernels do not each invent their
own constants or index mappings.

## Selected References

- Bisbas et al. 2012, original 3D-PDR architecture and model assumptions,
  MNRAS 427, 2100, https://doi.org/10.1111/j.1365-2966.2012.22077.x
- Tielens & Hollenbach 1985, common PDR notation for density, FUV field,
  chemistry, and thermal balance, ApJ 291, 722,
  https://doi.org/10.1086/163111
