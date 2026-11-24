---
title: "Hierarchy of Plane-Stress Elliptical Cavity Boundary Conditions"
domain: complex-elasticity
source_book: "Muskhelishvili / Inglis / Elliptical Cavity Progression"
section: "9-Level BC Waterfall"
tags:
  - muskhelishvili
  - inglis
  - boundary-conditions
  - julia
  - sympy
  - matlab
archive_source: "md/Gemini-Elliptical Cavity with Boundary displacements.md"
status: synthesized
---

# Hierarchy of Plane-Stress Elliptical Cavity Boundary Conditions

A structured progression of 2D plane-stress elliptical cavity problems ordered from simplest traction-free holes to general Fourier-prescribed displacements, designed around a shared `Common.jl` solver/plotting core:

1. **Traction-Free Elliptical Hole** (Inglis 1913 / Muskhelishvili): Zero inner wall tractions under uniaxial/biaxial far-field stress.
2. **Uniform Internal Pressure**: Constant normal pressure $P_{in}$ on cavity wall (superposed with far-field stress).
3. **Uniform Tangential Shear**: Constant tangential shear along the cavity rim.
4. **Rigid Inclusion (Perfect Bonding)**: $u = 0$ (or rigid-body rotation/translation) at cavity interface.
5. **Mixed Boundary Conditions**: Tractions on part of the rim, displacements on the rest (e.g., pin-loaded hole contact).
6. **Constant Normal Displacement (Case 1 Special)**: Uniform normal wall expansion $\rho = \text{const}$.
7. **Variable Normal Displacement (Case 1 General)**: Deforming from initial ellipse $(a_0, b_0)$ to $(a_1, b_1)$ via $\rho(\vartheta)$.
8. **Radial Displacement (Case 2)**: Outward radial displacement $l(\theta)$ from the ellipse center (dilatometer-like loading).
9. **Arbitrary Displacement Field**: Full Fourier series boundary displacement $u_x + i u_y = g(\varphi)$ (Eq. 34 in [[Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities]]).

## Symbolic Verification & Gotchas (Kirsch / Inglis Sanity Checks)

- **Exact vs. Numerical Derivatives**: Finite-differencing $\Phi(\zeta)$ near the unit circle $|\zeta|=1$ causes spurious non-zero normal/shear tractions ($T_x \neq 0$) on traction-free boundaries; always differentiate $\phi'(\zeta), \phi''(\zeta), \psi'(\zeta)$ analytically.
- **MATLAB `subs` Cell Array Gotcha**: Passing `{R, m, Gamma, GammaP, zeta}` vs `[...]` in Symbolic Math Toolbox requires matching flat vectors/structs when substituting complex parameters.
- **Multi-Language Kirsch Deduction Chain**: Ported across MATLAB Symbolic Math Toolbox, Python `SymPy` (`# %%` cells), and Julia `Symbolics.jl` confirming circular hole peak hoop stress $\sigma_{\theta\theta}(R, \pi/2) = 3\sigma_\infty$.

## Connections

- [[Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities]]
- [[Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole]]
- [[Hole Shape Optimization and Inverse Problems in MATLAB]]
