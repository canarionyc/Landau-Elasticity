---
title: "Hole Shape Optimization and Inverse Problems in MATLAB"
domain: computational-mechanics
source_book: "MATLAB PDE & Optimization Toolboxes"
section: "Forward vs Inverse Elasticity"
tags:
  - matlab
  - pde-toolbox
  - optimization
  - fmincon
  - inverse-problem
  - shape-optimization
archive_source: "md/Gemini-Elliptical Cavity with Boundary displacements.md"
status: synthesized
---

# Hole Shape Optimization and Inverse Problems in MATLAB

## Three MATLAB Approaches for the Kirsch Forward Problem

1. **Symbolic Math Toolbox**: Exact closed-form Kirsch polar stresses ($\sigma_{rr}, \sigma_{\theta\theta}, \tau_{r\theta}$) verifying $\text{SCF} = 3.0$.
2. **PDE Toolbox — Structural Mechanics (`createpde('structural','static-planestress')`)**: Industry-standard 2D plane-stress FEM (`structuralProperties`, `structuralBC`, `structuralLoads`, `res.VonMisesStress`).
3. **PDE Toolbox — General Coefficient Form (`createpde(2)`)**: Raw Navier-Cauchy PDE via the $c$-coefficient matrix of $(\lambda, G)$, useful for coupled multiphysics (thermoelasticity, non-standard laws).

## Inverse Problem: Shape Optimization with `fmincon`

Coupling `fmincon` with the PDE Toolbox minimizes peak Von Mises stress over ellipse semi-axes $(a, b)$ under finite-plate, multi-hole, or anisotropic conditions where analytical formulas do not exist. Under uniaxial $X$-tension, analytical hoop peak $\sigma_{\max} = \sigma_\infty(1 + 2b/a)$ drives the optimizer to elongate along $X$ ($b/a \to 0$).

## Preventing Singularities & Unphysical Boundary Pinning

- **Aspect Ratio Nonlinear Constraint**: Enforce $a/b \le 10$ (or $5$) via `nonlcon` so $b \not\to 0$ (which collapses the ellipse into a crack singularity).
- **Manufacturing ("Drill Bit") Floor**: Lower bound $R_{\min} \ge r_{\text{tool}}$.
- **Edge-Distance Keep-Out**: $L_{\text{edge}} \ge 2 \times \text{Hole Diameter}$ to prevent edge breakout.
- **Target Area Constraint**: $\pi a b = A_{\text{target}}$ prevents the optimizer from "cheating" by shrinking the hole to zero area.
- **Soft Barrier Penalty**: Add $1/(b - b_{\min})$ if design variables pin against box bounds `lb`.

## Connections

- [[Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole]]
- [[Hierarchy of Plane-Stress Elliptical Cavity Boundary Conditions]]
- [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]
