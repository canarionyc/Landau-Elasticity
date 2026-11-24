---
title: "Shell Edge Bending Boundary Layer (Geckeler Approximation)"
domain: shells-and-architecture
source_book: "Landau & Lifshitz / Shell Theory"
section: "§14 (Edge Bending)"
tags:
  - landau
  - shells
  - boundary-layer
  - geckeler-approximation
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Shell Edge Bending Boundary Layer (Geckeler Approximation)

When an edge ring or tension chain restrains the outward radial expansion predicted by [[Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)|pure membrane theory]], bending stiffness $D$ must be restored locally near the base ($y = 0$, measuring meridional distance $y$ upward along the shell).

## Governing ODE: Beam on an Elastic Foundation

Outward normal deflection $w(y)$ is resisted by meridional wall bending $D$ and circumferential hoop membrane stretching $k$:

$$
D \frac{d^4 w}{dy^4} + k w = 0, \qquad D = \frac{Eh^3}{12(1-\nu^2)}, \qquad k = \frac{Eh}{a^2}
$$

## The Great Material Cancellation

Defining $4\beta^4 = k/D$, Young's modulus $E$ cancels identically between bending and hoop stiffness:

$$
4\beta^4 = \frac{Eh / a^2}{Eh^3 / [12(1-\nu^2)]} = \frac{12(1-\nu^2)}{a^2 h^2} \implies \beta = \frac{\left[3(1-\nu^2)\right]^{1/4}}{\sqrt{a h}}
$$

Thus, the boundary layer decay length $1/\beta \sim \sqrt{ah}$ is **purely geometric**, independent of $E$!

## Decaying Oscillatory Solution

$$
w(y) = e^{-\beta y}\left(C_1 \cos(\beta y) + C_2 \sin(\beta y)\right)
$$

Bending moments induced by the edge restraint decay exponentially away from the base, recovering pure membrane behavior further up the dome.

## Connections

- [[Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
- [[Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)]]
