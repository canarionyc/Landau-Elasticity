---
title: "FDM vs FEM and Hermite Cubic Shape Functions for Buckling"
domain: numerical-methods
source_book: "Landau & Lifshitz §21 / Computational Mechanics"
section: "Strong Form vs Weak Form"
tags:
  - fdm
  - fem
  - buckling
  - hermite-polynomials
  - julia
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/Clamped-Clamped Buckling (FDM vs FEM).jl"
status: synthesized
---

# FDM vs FEM and Hermite Cubic Shape Functions for Buckling

Comparing **Finite Difference (FDM)** and **Finite Element (FEM)** on the clamped-clamped 1D Euler buckling problem $EI w'''' + P w'' = 0$ ($P_{cr} = \frac{4\pi^2 EI}{L^2}$):

| Aspect | Finite Difference (FDM) | Finite Element (FEM) |
| --- | --- | --- |
| **Approximates** | Differential operators (Taylor stencils) | Solution space (piecewise basis functions) |
| **Formulation** | Strong form ($w'''', w''$ directly) | Weak form via integration by parts: $\int EI w'' v''\,dx = P\int w' v'\,dx$ |
| **Max Derivative** | 4th derivative (`[1,-4,6,-4,1]`) | 2nd derivative ($C^1$ cubic Hermite shape functions) |
| **Accuracy ($N=6$)** | ~5–10% error | $< 0.5\%$ error |
| **Best For** | Simple rectangular grids | Complex 3D frames, arbitrary geometries, mixed BCs |

See [Clamped-Clamped Buckling (FDM vs FEM).jl](../../examples/Clamped-Clamped%20Buckling%20%28FDM%20vs%20FEM%29.jl).

## Subtlety: Why Cubic FEM is Approximate for Buckling but Exact for Static Point Loads

- **Buckling ($EI w'''' + |T|w'' = 0$)**: Exact solution is **transcendental/trigonometric** ($A + Bx + C\sin(kx) + D\cos(kx)$). Piecewise cubic Hermite elements only approximate the sine wave, converging rapidly as mesh refines.
- **Static Bending under Point Loads ($EI w'''' = 0$ between loads)**: Bending moment $M(x)$ is linear $\implies w''$ is linear $\implies w(x)$ is **piecewise cubic**. Here, cubic Hermite FEM reproduces the exact analytical deflection with minimal elements.

## Connections

- [[Small Rod Deflection and Euler Buckling (Landau §20-§21)]]
- [[Plate Buckling and Garland Curves in Box Girder Bridges]]
- [[Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev]]
- [[2D Truss and Frame Data Structures in Rust]]
- [[3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly]]
