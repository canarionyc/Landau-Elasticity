---
title: "Torsion of Non-Circular Rods and Prandtl Soap-Film Analogy (Landau §16)"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§16"
tags:
  - landau
  - torsion
  - prandtl-stress-function
  - warping
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/Ex16.2 - Torsion of an Elliptical Rod (§16, Problem 2).jl"
  - "examples/Ex16.3 - Torsion of a Triangular Rod (§16, Problem 3).jl"
status: synthesized
---

# Torsion of Non-Circular Rods and Prandtl Soap-Film Analogy (Landau §16)

Unlike circular shafts where plane cross-sections remain flat, non-circular cross-sections **warp** out of plane under torsion: $u_z = \tau \psi(x, y)$.

## Prandtl Stress Function

Introducing the stress function $\chi(x, y)$ (such that shear stresses are $\sigma_{xz} = \partial\chi/\partial y$ and $\sigma_{yz} = -\partial\chi/\partial x$) reduces torsion to a Poisson problem with Dirichlet boundary condition $\chi = 0$ on the cross-section perimeter:

$$
\nabla^2 \chi = -2\mu\tau
$$

The torsional rigidity is twice the integral of $\chi$ over the cross-section: $C = 2 \iint \chi\, dx\, dy$.

## Exact Cross-Section Examples

- **Elliptical Rod** (*§16 Problem 2*): $\chi(x,y) \propto (1 - x^2/a^2 - y^2/b^2)$, giving $C = \frac{\pi \mu a^3 b^3}{a^2 + b^2}$. See [Ex16.2 - Torsion of an Elliptical Rod (§16, Problem 2).jl](../../examples/Ex16.2%20-%20Torsion%20of%20an%20Elliptical%20Rod%20%28%C2%A716,%20Problem%202%29.jl).
- **Equilateral Triangular Rod** (*§16 Problem 3*): Cubic polynomial $\chi$ vanishing on all three edges. See [Ex16.3 - Torsion of a Triangular Rod (§16, Problem 3).jl](../../examples/Ex16.3%20-%20Torsion%20of%20a%20Triangular%20Rod%20%28%C2%A716,%20Problem%203%29.jl).

## Prandtl's Soap-Film Analogy & Why I-Beams Twist Easily

Inflating a membrane (soap film) stretched over a hole shaped like the cross-section obeys the exact same Poisson equation:
- **Bubble slope** $\leftrightarrow$ local shear stress magnitude.
- **Volume under the bubble** $\leftrightarrow$ torsional rigidity $C$ (or engineering $J$).

For thin open sections like an **I-beam** (narrow rectangular web and flanges), the soap film barely rises, so $J \approx \frac{1}{3}\sum b t^3$ is tiny compared to $I_z$ or $I_y$ (see [[3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly]]). Closed hollow tubes trap a large bubble volume and are orders of magnitude stiffer in torsion.

## Connections

- [[Helicoidal Rods, Coil Springs, and Elastic Free Energy (Landau §19)]]
- [[3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly]]
- [[Modernizing Love's Elasticity with 4th-Order Tensors]]
- [[Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev]] (Poisson equation $\nabla^2 u = f$ with Dirichlet BCs)
