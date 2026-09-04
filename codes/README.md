# Codes

MATLAB scripts for the vibration analysis of spacecraft.

Lumped mass models.

## [vibroacoustics.m](vibroacoustics.m)

Chap. 2, [slides/svs_02_dinamica_1gdl.pdf](../slides/svs_02_dinamica_1gdl.pdf).

Single-dof plate model with a dummy mass subjected to acoustic loadings. Given the launcher SPL, it computes:
- pressure spectral density;
- transfer function;
- PSD of the dummy mass acceleration;
- PSD of the displacement;
- rms values.

## [launcher.m](launcher.m)

Chap. 3, [slides/svs_03_dinamica_mgdl.pdf](../slides/svs_03_dinamica_mgdl.pdf).

For a 4-dof payload-launcher system, computes:
- mode shapes and natural frequencies;
- dynamic response analysis;
- effective modal masses and modal participation factors;
- effect of each mode on the response.

## [guyan.m](guyan.m)

Chap. 4, [slides/svs_04_dinamica_substructuring.pdf](../slides/svs_04_dinamica_substructuring.pdf).

For a 10-dof payload-launcher system, computes:
- mode shapes and natural frequencies;
- generalized mass and stiffness matrices;
- statically condensation (Guyan reduction);
- reduced mass and stiffness;
- MAC and COC.

## [effmasses.m](effmasses.m)*

Additional example Chap. 4, credits to [@PieroChiaia](https://github.com/PieroChiaia)

For the same 10-dof model of the previous example, computes:
- mode shape animations;
- generalized mass and stiffness matrices;
- participation factors;
- modal effective masses.

*Function [animate.m](animate.m) is needed to run the script.

## [buckling.m](buckling.m)

Chap. 10-12, [slides/](../slides/). 

For a column in compression, computes:
- the first three critical buckling loads;
- the first three buckling mode shapes.

## [postbuckling.m](postbuckling.m)

For the same problem of above, computes:
- the equilibrium curve in the post-buckling state;
- comparison with the linearized solution.

## [elastica.nb](elastica.nb)

Mathematica notebook. Same exercise as in [postbuckling.m](postbuckling.m).

## [pendulum.nb](pendulum.nb)

Non linear pendulum; analogy with the elastica solution.

## [solararray.m](solararray.m)

Chap. 12, [slides/svs_12_esempi_nonlineari.pdf](../slides/svs_12_esempi_nonlineari.pdf).

Kinematic analysis of a multi-body system representing a deployable solar array.
It provides:
- animation states of the deployed structure
- actuation mechanics

##[tapedeploy.nb](tapedeploy.nb)

Chap. 12, [slides/svs_12_esempi_nonlineari.pdf](../slides/svs_12_esempi_nonlineari.pdf).

Single dof deployment analysis of a tape spring. Wolfram Mathematica notebook.
It provides:
- full formulation of the unconstrained expanding circle model
- dissipation and gravity effects
- deployment speed and dynamics
- comparison with experiment
