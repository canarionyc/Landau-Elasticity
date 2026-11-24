---
title: "Small Rod Deflection and Euler Buckling (Landau §20-§21)"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§20, §21"
tags:
  - landau
  - buckling
  - thin-rods
  - stability
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/Small Deflections of a Clamped-Supported Rod (§20).jl"
  - "examples/The Stability of Elastic Systems - Euler Buckling (§21).jl"
  - "examples/The Stability of Elastic Systems - Clamped Ends (§21).jl"
status: synthesized
---

# Small Rod Deflection and Euler Buckling (Landau §20-§21)

## Small Deflections (§20)

For small transverse displacements $w(x)$ under transverse load $q(x)$, Euler-Bernoulli beam equilibrium gives:

$$
EI \frac{d^4 w}{dx^4} = q(x)
$$

See [Small Deflections of a Clamped-Supported Rod (§20).jl](../../examples/Small%20Deflections%20of%20a%20Clamped-Supported%20Rod%20%28%C2%A720%29.jl).

## Euler Buckling (§21)

Under an axial compressive load $T$, linearized stability is governed by:

$$
EI \frac{d^4 w}{dx^4} + T \frac{d^2 w}{dx^2} = 0
$$

The general solution is trigonometric + linear: $w(x) = A + Bx + C\sin(kx) + D\cos(kx)$ with $k = \sqrt{T/EI}$ (see [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]).

- **Hinged (Pinned-Pinned) Ends**: $T_{cr} = \frac{\pi^2 EI}{L^2}$. See [The Stability of Elastic Systems - Euler Buckling (§21).jl](../../examples/The%20Stability%20of%20Elastic%20Systems%20-%20Euler%20Buckling%20%28%C2%A721%29.jl).
- **Clamped-Clamped Ends**: Restraining end rotations doubles the effective wavenumber and increases critical buckling load by **4x**:

$$
T_{clamped} = \frac{4\pi^2 EI}{L^2}
$$

See [The Stability of Elastic Systems - Clamped Ends (§21).jl](../../examples/The%20Stability%20of%20Elastic%20Systems%20-%20Clamped%20Ends%20%28%C2%A721%29.jl).

## Limitation of Linear Buckling

Linear theory leaves post-buckling amplitude indeterminate at $T = T_{cr}$ (neutral horizontal line on load-deflection plot). Resolving post-buckling requires [[Euler Elastica and Post-Buckling Stability]].

## Connections

- [[Plate Buckling and Garland Curves in Box Girder Bridges]]
- [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]
- [[Euler Elastica and Post-Buckling Stability]]
- [[2D Truss and Frame Data Structures in Rust]]
