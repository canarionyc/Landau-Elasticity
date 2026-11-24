---
title: "Helicoidal Rods, Coil Springs, and Elastic Free Energy (Landau §19)"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§19"
tags:
  - landau
  - thin-rods
  - springs
  - free-energy
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/Ex19.6 - Rod Twisted into a Spiral (§19, Problem 6).jl"
status: synthesized
---

# Helicoidal Rods, Coil Springs, and Elastic Free Energy (Landau §19)

In Kirchhoff thin-rod theory, a straight inextensible wire of length $L$ wound into a helix of radius $R$ and pitch angle $\alpha$ maintains material length $L$ while its macroscopic height shortens geometrically to:

$$
H = L \sin\alpha, \qquad N = \frac{L\cos\alpha}{2\pi R}
$$

Maintaining equilibrium requires a purely axial force $F_z$ **plus** balancing twisting and bending moments ($M_z, M_\phi$) at the ends. See [Ex19.6 - Rod Twisted into a Spiral (§19, Problem 6).jl](../../examples/Ex19.6%20-%20Rod%20Twisted%20into%20a%20Spiral%20%28%C2%A719,%20Problem%206%29.jl).

## Elastic Free Energy vs. Classical Spring Potential Energy

For an isothermal, reversible deformation, classical mechanical potential energy $U$ equals the thermodynamic **Helmholtz Free Energy** $F - F_0$:

$$
U = \int_0^L \left(\frac{1}{2} E I_x \kappa_x^2 + \frac{1}{2} E I_y \kappa_y^2 + \frac{1}{2} C \tau^2\right) dl
$$

For a standard closely-coiled helical spring (small pitch $\alpha$), axial deflection $\Delta z$ is accommodated almost entirely by **torsion** of the wire (~90–95% of stored energy), yielding $U = \frac{1}{2} k (\Delta z)^2$ with spring constant governed by the shear modulus $\mu$ (not $E$):

$$
k = \frac{C}{2\pi N R^3} = \frac{\mu r^4}{4 N R^3} = \frac{G d^4}{64 N R^3}, \qquad C = \frac{1}{2}\mu\pi r^4, \qquad \mu = \frac{E}{2(1+\sigma)}
$$

## Engineering & Biological Applications

- **Mechanical Coil Springs**: Convert stiff metal shear rigidity into large axial compliance.
- **Helical Antennas**: Chiral geometry radiates circularly polarized waves resistant to ionospheric Faraday rotation.
- **Bacterial Flagella & Micro-robots**: Convert rotation in viscous fluid into axial thrust.
- **Wire Ropes & Heat Exchangers**: Flexibility via strand slip; Dean vortices in coiled pipes enhance heat transfer.

## Connections

- [[Torsion of Non-Circular Rods and Prandtl Soap-Film Analogy (Landau §16)]]
- [[Small Rod Deflection and Euler Buckling (Landau §20-§21)]]
- [[Bifurcation Theory and Double-Well Free Energy]]
