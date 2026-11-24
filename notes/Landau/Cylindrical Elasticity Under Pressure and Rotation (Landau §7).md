---
title: "Cylindrical Elasticity Under Pressure and Rotation (Landau §7)"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§7"
tags:
  - landau
  - elasticity
  - cylindrical-symmetry
  - pressure-vessel
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/Example 4 - Cylindrical Pipe Under Internal Pressure.jl"
  - "examples/Example 5 - Rotating Cylinder.jl"
status: synthesized
---

# Cylindrical Elasticity Under Pressure and Rotation (Landau §7)

In cylindrical coordinates $(r, \phi, z)$ with radial symmetry ($u_r = u(r)$, $u_\phi = 0$), the homogeneous equilibrium equation gives:

$$
u(r) = A r + \frac{B}{r}
$$

## Canonical Problems

1. **Thick-Walled Cylindrical Pipe Under Internal Pressure** (*Example 4*): With $\sigma_{rr}(R_1) = -P$ and $\sigma_{rr}(R_2) = 0$, the radial stress is compressive while the hoop stress $\sigma_{\phi\phi}$ is tensile and peaks at the inner wall $r = R_1$. See [Example 4 - Cylindrical Pipe Under Internal Pressure.jl](../../examples/Example%204%20-%20Cylindrical%20Pipe%20Under%20Internal%20Pressure.jl).
2. **Rotating Solid Cylinder** (*Example 5*): Centrifugal body force $\rho \omega^2 r$ adds a particular solution proportional to $r^3$, bounded at $r = 0$ and traction-free at $r = R$. See [Example 5 - Rotating Cylinder.jl](../../examples/Example%205%20-%20Rotating%20Cylinder.jl).

## Thin-Wall vs. Thick-Wall Limit

As $R_2 - R_1 = h \ll R$, the radial stress varies negligibly compared to the hoop stress, recovering Mariotte's membrane hoop tension $N_\phi = P R$ (see [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]). For thick walls, inner and outer pressures must remain separate boundary conditions rather than a single net resultant.

## Connections

- [[Spherically Symmetric Elasticity (Landau §7)]]
- [[Free-Form Architecture vs Symmetric Coordinates (Morpheus Case)]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
