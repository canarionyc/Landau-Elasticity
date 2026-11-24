---
title: "Plate Buckling and Garland Curves in Box Girder Bridges"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§21 (2D Plate Extension)"
tags:
  - landau
  - plate-buckling
  - bridges
  - fdm
  - garland-curves
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
related_code:
  - "examples/The Garland Curve.jl"
  - "examples/Clamped Plate Buckling.jl"
status: synthesized
---

# Plate Buckling and Garland Curves in Box Girder Bridges

Extending 1D rod buckling to a 2D compressed plate (such as the top flange of a steel box girder bridge) with flexural rigidity $D = \frac{Eh^3}{12(1-\sigma^2)}$ gives:

$$
D\left(\frac{\partial^4 \zeta}{\partial x^4} + 2\frac{\partial^4 \zeta}{\partial x^2\partial y^2} + \frac{\partial^4 \zeta}{\partial y^4}\right) + N_x \frac{\partial^2 \zeta}{\partial x^2} = 0
$$

## Simply Supported Analytical Solution (Garland Curves)

For a simply supported plate of length $a$ and width $b$ buckling into $m$ half-waves along $x$:

$$
\sigma_{cr} = \frac{k \pi^2 E}{12(1-\sigma^2)}\left(\frac{h}{b}\right)^2, \qquad k = \left(\frac{mb}{a} + \frac{a}{mb}\right)^2
$$

Plotting $k(a/b)$ for $m = 1, 2, 3, \dots$ produces the **Garland Curves** whose lower envelope never drops below $k_{\min} = 4.0$ for long plates ($a \gg b$). Welding longitudinal **stiffeners** reduces effective panel width $b$, boosting $\sigma_{cr} \propto (h/b)^2$ above steel yield strength. See [The Garland Curve.jl](../../examples/The%20Garland%20Curve.jl).

## Clamped Plate Buckling via Finite Differences (FDM)

For clamped edges ($\zeta = 0, \partial\zeta/\partial n = 0$), Fourier modes couple and analytical solutions become cumbersome. Discretizing with 1D stencils ($D_{2x}, D_{4x}$ with ghost-point modification $u_{-1} = u_1 \implies$ diagonal $+1$) and Kronecker products turns the biharmonic operator into:

$$
L_{\text{biharm}} = D_{4x} \otimes I_y + 2(D_{2x} \otimes D_{2y}) + I_x \otimes D_{4y}
$$

Solving the generalized eigenvalue problem $L_{\text{biharm}} u = \lambda (-D_{2x} \otimes I_y) u$ raises $k$ from $4.0$ (simply supported) to $\approx 6.9\text{--}7.2$ at $a/b = 2$. See [Clamped Plate Buckling.jl](../../examples/Clamped%20Plate%20Buckling.jl).

## Connections

- [[Small Rod Deflection and Euler Buckling (Landau §20-§21)]]
- [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]
- [[Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)]]
