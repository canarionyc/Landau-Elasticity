---
title: "Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole"
domain: complex-elasticity
source_book: "Kirsch / Muskhelishvili / Mansfield (1953)"
section: "Optimal & Neutral Holes"
tags:
  - kirsch
  - muskhelishvili
  - mansfield
  - neutral-hole
  - stress-concentration
  - matlab
archive_source: "md/Gemini-Muskhelishvili Elasticity_ Complex Variables.md"
status: synthesized
---

# Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole

Under biaxial far-field tension $(S_x, S_y)$ (e.g. $2:1$ fuselage ratio $S_x = 200\text{ MPa}, S_y = 100\text{ MPa}$), cutout stress control progresses through three classical stages:

1. **Unreinforced Circular Hole (Kirsch)**:
   Hoop stress around the rim varies as $\sigma_\theta(\theta) = (S_x + S_y) - 2(S_x - S_y)\cos(2\theta)$. Produces a "butterfly" stress field with sharp red hotspots ($\sigma > 2.5\times$ load; $\text{SCF}=3$ under uniaxial load) and under-stressed zones.
2. **Harmonic / Optimal Ellipse (Muskhelishvili Inverse Solution)**:
   Choosing an ellipse with semi-axes matched to the far-field stress ratio $a/b = S_x/S_y$ eliminates stress peaks along the boundary, achieving **constant hoop stress** around the entire rim.
3. **Mansfield Neutral ("Invisible") Hole**:
   Adding a variable-thickness reinforcement liner along the ellipse rim balances the traction normal to the cut ($T_n = S_x\sin^2\theta + S_y\cos^2\theta$), leaving the surrounding plate **completely undisturbed** with uniform stress everywhere outside the liner.

## MATLAB Visualization Workflow (2x2 Layout)

`Gemini-Muskhelishvili Elasticity_ Complex Variables.md` contains the full self-contained MATLAB script plotting:
- Top-left: Von Mises contour for Kirsch circular hole.
- Top-right: Harmonic ellipse ($a/b = S_x/S_y$) with constant boundary hoop stress.
- Bottom-left: Mansfield neutral hole with variable-thickness black reinforcement liner.
- Bottom-right: Polar comparison of boundary $\sigma_\theta(\theta)$ and Mansfield liner thickness.

## Connections

- [[Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities]]
- [[Hole Shape Optimization and Inverse Problems in MATLAB]]
- [[Spherically Symmetric Elasticity (Landau §7)]] (3D spherical cavity analogue)
- [[Viscoelastic Stress Relaxation - Zener vs Rabotnov Models (Savin Ch. XII)]]
