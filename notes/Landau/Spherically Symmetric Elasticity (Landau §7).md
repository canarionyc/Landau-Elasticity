---
title: "Spherically Symmetric Elasticity (Landau §7)"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§7"
tags:
  - landau
  - elasticity
  - spherical-symmetry
  - stress-concentration
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/Example 3 - Self-Gravitating Sphere.jl"
  - "examples/Example 12 - Stress Concentration around a Spherical Cavity.jl"
status: synthesized
---

# Spherically Symmetric Elasticity (Landau §7)

For purely radial displacement $\mathbf{u} = u(r)\hat{\mathbf{r}}$ in an isotropic body with zero curl ($\nabla \times \mathbf{u} = 0$), the equilibrium equation reduces to:

$$
\frac{d}{dr}\left(\frac{1}{r^2}\frac{d(r^2 u)}{dr}\right) = -\frac{\rho f_r}{\lambda + 2\mu}
$$

Without body forces, the general radial displacement is:

$$
u(r) = A r + \frac{B}{r^2}
$$

## Canonical §7 Problems

1. **Solid Sphere in Own Gravitational Field** (*Example 3*): Regularity at $r = 0$ removes the $B/r^2$ singularity. Surface traction vanishes at $r = R$, yielding central compression and tension/compression transition controlled by Poisson's ratio $\sigma$. See [Example 3 - Self-Gravitating Sphere.jl](../../examples/Example%203%20-%20Self-Gravitating%20Sphere.jl).
2. **Hollow Sphere Under Internal and External Pressure** (*§7 Problem 1*): Boundary conditions $\sigma_{rr}(a) = -p_1$ and $\sigma_{rr}(b) = -p_2$ determine both $A$ and $B$ (Lamé's thick-walled sphere). Contrast this with [[Free-Form Architecture vs Symmetric Coordinates (Morpheus Case)|thick-solid vs. thin-shell pressure loading]].
3. **Spherical Cavity in an Infinite Medium** (*§7 Problem 2 / Example 12*): Under hydrostatic tension at infinity or internal cavity pressure, $A$ and $B$ produce local stress concentration around a void. See [Example 12 - Stress Concentration around a Spherical Cavity.jl](../../examples/Example%2012%20-%20Stress%20Concentration%20around%20a%20Spherical%20Cavity.jl).

## Modern Tensor Verification

Using automated symbolic calculus via [[From Love Scale Factors to Metric Tensors and Christoffel Symbols]], substituting the spherical metric $g_{ij} = \mathrm{diag}(1, r^2, r^2\sin^2\theta)$ recovers the exact Navier radial equilibrium equation:

$$
(\lambda + 2\mu)\left(\frac{d^2 u}{dr^2} + \frac{2}{r}\frac{du}{dr} - \frac{2u}{r^2}\right) = 0
$$

## Connections

- [[Cylindrical Elasticity Under Pressure and Rotation (Landau §7)]]
- [[From Love Scale Factors to Metric Tensors and Christoffel Symbols]]
- [[Finite Elasticity and Anisotropic Curved Solids]] (nonlinear Neo-Hookean spherical shell inflation)
- [[Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole]] (2D cavity analogue)
