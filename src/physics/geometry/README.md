# Geometry

This directory contains HEALPix utilities, ray-path helpers, and evaluation
point construction. It turns point-cloud geometry into the per-ray distances,
projected point ids, and adaptive stepping information needed by column,
radiation, and LVG calculations.

## Background

3D-PDR estimates directional attenuation by casting rays through the modeled
cloud. HEALPix provides an equal-area angular sampling of the sphere, which is
useful when averaging radiation or escape probabilities over many directions.
Evaluation points then sample material along each ray so columns, visual
extinction, and optical depths can be integrated.

## Selected References

- Gorski et al. 2005, HEALPix angular pixelization, ApJ 622, 759,
  https://doi.org/10.1086/427976
- Bisbas et al. 2012, 3D-PDR ray tracing and three-dimensional PDR treatment,
  MNRAS 427, 2100, https://doi.org/10.1111/j.1365-2966.2012.22077.x
