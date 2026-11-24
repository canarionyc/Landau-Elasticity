---
title: "Toroidal Coordinates and the Grad-Shafranov Equation"
domain: tensor-calculus-and-plasma
source_book: "Curvilinear Tensor Calculus / Tokamak Equilibrium"
section: "Toroidal Metric"
tags:
  - toroidal-coordinates
  - christoffel-symbols
  - grad-shafranov
  - fusion
  - laplacian
archive_source: "md/Gemini-AEH Love Elasticity Theory_ Classical to Modern.md"
status: synthesized
---

# Toroidal Coordinates and the Grad-Shafranov Equation

Swapping the spherical metric in [[From Love Scale Factors to Metric Tensors and Christoffel Symbols]] for a **toroidal metric** $(r, \theta, \zeta)$ models Tokamak fusion reactors (ITER, JET) where magnetic pressure balances $\nabla p = \mathbf{J} \times \mathbf{B}$.

## Toroidal Metric Tensor

With major radius $R_0$, minor radius $r$, poloidal angle $\theta$, and toroidal angle $\zeta$:

$$
h_r = 1, \qquad h_\theta = r, \qquad h_\zeta = R_0 + r\cos\theta \implies ds^2 = dr^2 + r^2 d\theta^2 + (R_0 + r\cos\theta)^2 d\zeta^2
$$

Computing $\nabla^2 P = \frac{1}{\sqrt{g}}\partial_i(\sqrt{g}\, g^{ij}\partial_j P)$ under toroidal symmetry ($\partial_\zeta P = 0$) yields cylindrical terms plus **toroidal coupling terms** proportional to $\frac{1}{R_0 + r\cos\theta}$:

$$
\nabla^2 P = \frac{\partial^2 P}{\partial r^2} + \frac{1}{r}\frac{\partial P}{\partial r} + \frac{1}{r^2}\frac{\partial^2 P}{\partial\theta^2} + \frac{\cos\theta}{R_0 + r\cos\theta}\frac{\partial P}{\partial r} - \frac{\sin\theta}{r(R_0 + r\cos\theta)}\frac{\partial P}{\partial\theta}
$$

## Physical Consequence: Shafranov Shift & Grad-Shafranov

Because the inner side of the torus is tighter than the outer side, the Christoffel term $\Gamma^r_{\zeta\zeta}$ produces an outward radial force pushing plasma flux surfaces toward $\theta = 0$ (**Shafranov Shift**). Equilibrium magnetic flux $\psi$ obeys the **Grad-Shafranov equation**:

$$
\Delta^* \psi = -\mu_0 R^2 \frac{dp}{d\psi} - F\frac{dF}{d\psi}
$$

requiring outer poloidal field coils and D-shaped cross-sections in modern reactors.

## Connections

- [[From Love Scale Factors to Metric Tensors and Christoffel Symbols]]
- [[Free-Form Architecture vs Symmetric Coordinates (Morpheus Case)]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
