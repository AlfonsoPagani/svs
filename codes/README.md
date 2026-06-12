# Codes

MATLAB scripts for the vibration analysis of spacecraft.

Lumped mass models.

## [vibroacoustics.m](vibroacoustics.m)

Sect. 2.9, [Lecture notes](book/svsbook.pdf).

Single-dof plate model with a dummy mass subjected to acoustic loadings. Given the launcher SPL, it computes:
- pressure spectral density;
- transfer function;
- PSD of the dummy mass acceleration;
- PSD of the displacement;
- rms values.

## [launcher.m](launcher.m)

Chap. 3, [Lecture notes](book/svsbook.pdf).

For a 4-dof payload-launcher system, computes:
- mode shapes and natural frequencies;
- dynamic response analysis;
- effective modal masses and modal participation factors;
- effect of each mode on the response.

## [guyan.m](guyan.m)

Sect. 4.1, [Lecture notes](book/svsbook.pdf).

For a 10-dof payload-launcher system, computes:
- mode shapes and natural frequencies;
- generalized mass and stiffness matrices;
- statically condensation (Guyan reduction);
- reduced mass and stiffness;
- MAC and COC.
