---
title: "General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)"
domain: shells-and-tensors
source_book: "A. E. H. Love, Treatise on the Mathematical Theory of Elasticity"
section: "Ch. XXI-XXIV (Arts. 322, 331)"
tags:
  - love
  - shells
  - fundamental-forms
  - paraboloid
  - curvature-tensor
archive_source: "md/Gemini-AEH Love Elasticity Theory_ Classical to Modern.md"
status: synthesized
---

# General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)

A thin shell is a 2D surface $\mathbf{r}(u, v)$ embedded in 3D space, governed by **two** surface tensors:
1. **First Fundamental Form (Metric $g_{ij}$)**: Measures in-surface stretching ($g_{ij} = \mathbf{e}_i \cdot \mathbf{e}_j$).
2. **Second Fundamental Form (Curvature $b_{ij}$)**: Measures out-of-plane bending via the unit normal $\mathbf{n}$ ($b_{ij} = \mathbf{n} \cdot \partial_j \mathbf{e}_i$).

## Axisymmetric Paraboloid Membrane Solution

For a paraboloid of revolution $z = u^2/a$ under uniform pressure $P$, axisymmetry sets $\partial/\partial v = 0$ and $T_{uv} = 0$, reducing the PDEs to algebraic/ODE balance using principal radii of curvature $R_1 = g_{11}/b_{11}$ and $R_2 = g_{22}/b_{22}$:

1. **Normal Equilibrium (Laplace / Soap-Bubble Equation)**:

$$
\frac{T_1}{R_1} + \frac{T_2}{R_2} = P
$$

2. **Global Cap Balance**: $T_1 (2\pi u)\sin\phi = P(\pi u^2)$, yielding exact analytical meridional and hoop tensions:

$$
T_1 = \frac{Pa}{2}\sqrt{1 + \frac{4u^2}{a^2}}, \qquad T_2 = Pa\sqrt{1 + \frac{4u^2}{a^2}} - T_1
$$

## Geometric Subtlety: Dome Orientation & Hoop Curvature

- **Dome vs. Dish**: $z = u^2/a$ opens upward (dish); to model a dome under downward load or upward pressure, flip sign ($z = -u^2/a$) or reverse normal orientation.
- **Hoop Curvature $1/R_2 \neq 1/u$**: By Meusnier's theorem, the normal hoop radius of curvature $R_2$ is the distance along the surface normal to the symmetry axis, $R_2 = u / \sin\alpha$ (where $\tan\alpha = dz/du = 2u/a$). At the crown ($u \to 0$), $1/R_2 \to 2/a$ remains **finite**, whereas $1/u$ would falsely blow up.

## Analytical vs. Non-FEA Numerical Solution

For non-symmetric shells (e.g., elliptical domes), the tensor engine derives the **exact governing PDE** on the smooth surface, which can then be solved to exponential accuracy via Chebyshev spectral collocation (see [[Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev]]) without tetrahedral FEA meshing.

## Connections

- [[From Love Scale Factors to Metric Tensors and Christoffel Symbols]]
- [[Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)]]
- [[Shell Edge Bending Boundary Layer (Geckeler Approximation)]]
- [[Free-Form Architecture vs Symmetric Coordinates (Morpheus Case)]]
- [[Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev]]
