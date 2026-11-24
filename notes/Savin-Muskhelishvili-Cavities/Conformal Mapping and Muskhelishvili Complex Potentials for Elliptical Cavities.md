---
title: "Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities"
domain: complex-elasticity
source_book: "N. I. Muskhelishvili / Elliptical Cavity Paper"
section: "2D Complex Variable Method"
tags:
  - muskhelishvili
  - conformal-mapping
  - elliptical-cavity
  - complex-potentials
  - julia
archive_source: "md/Gemini-Elliptical Cavity with Boundary displacements.md"
status: synthesized
---

# Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities

In 2D plane elasticity, elliptical coordinates $(p, \vartheta)$ conveniently express boundary normal/tangential displacements, while **Muskhelishvili's conformal mapping** maps the exterior of an ellipse in the physical $z$-plane onto the exterior of the unit circle ($|\zeta| \ge 1$) in the $\zeta$-plane.

## 1. Conformal Mapping & Inverse

$$
z = \omega(\zeta) = R\left(\zeta + \frac{m}{\zeta}\right), \qquad R = \frac{a + b}{2}, \qquad m = \frac{a - b}{a + b}
$$

As $|z| \to \infty$, $z \approx R\zeta$, so far-field principal stresses $P_1, P_2$ (at angle $\Lambda$ to the $x$-axis) scale cleanly into potential asymptotics:

$$
\Gamma = -\frac{P_1 + P_2}{4}, \qquad \Gamma' = \frac{P_1 - P_2}{2}e^{-2i\Lambda}
$$

To invert $z \mapsto \zeta$ numerically, solve $R\zeta^2 - z\zeta + Rm = 0$ and select the root with $|\zeta| \ge 1$.

## 2. General Displacement-Controlled Potentials (Eq. 34)

For arbitrary boundary displacements expanded in Fourier coefficients $A_{\pm(2n-1)}$ (with shear modulus $G$ and Kolosov constant $\kappa$):

$$
\phi(\zeta) = \Gamma R\zeta + \frac{2G}{\kappa}\sum_{n=1}^{\infty} \frac{A_{-(2n-1)}}{\zeta^{2n-1}} + (m\Gamma + \Gamma')\frac{R}{\kappa\zeta}
$$

$$
\psi(\zeta) = \Gamma' R\zeta - 2G\sum_{n=1}^{\infty} \frac{A_{2n-1}}{\zeta^{2n-1}} + \Gamma R\left(\frac{\kappa}{\zeta} - \frac{\zeta(1+m^2)}{\zeta^2 - m}\right) - \frac{\zeta(1+m\zeta^2)}{\zeta^2 - m}\phi'(\zeta)_{\text{series}}
$$

## 3. Kolosov-Muskhelishvili Stress Formulas

Using $\Phi(\zeta) = \phi'(\zeta)/\omega'(\zeta)$ (evaluated via **exact analytical derivatives** rather than finite differences to avoid boundary traction noise):

$$
\sigma_x + \sigma_y = 4\,\mathrm{Re}\left[\frac{\phi'(\zeta)}{\omega'(\zeta)}\right]
$$

$$
\sigma_y - \sigma_x + 2i\tau_{xy} = 2\left[\frac{\overline{\omega(\zeta)}}{\omega'(\zeta)}\left(\frac{\phi'(\zeta)}{\omega'(\zeta)}\right)' + \frac{\psi'(\zeta)}{\omega'(\zeta)}\right]
$$

## Connections

- [[Hierarchy of Plane-Stress Elliptical Cavity Boundary Conditions]]
- [[Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole]]
- [[Hole Shape Optimization and Inverse Problems in MATLAB]]
- [[Modernizing Love's Elasticity with 4th-Order Tensors]]
