---
title: "Linear Stability Analysis and Jacobian Eigenvalues (Lotka-Volterra)"
domain: dynamical-systems
source_book: "Mathematical Ecology / Nonlinear Dynamics (Off-Topic)"
section: "Jacobian Linearization"
tags:
  - off-topic
  - dynamical-systems
  - jacobian
  - eigenvalues
  - stability
archive_source: "md/Gemini-Predator-Prey Models and Outcomes.md"
status: synthesized
---

# Linear Stability Analysis and Jacobian Eigenvalues (Lotka-Volterra)

*(Off-topic exploration outside elasticity, saved from Gemini chat archive.)*

Step-by-step linearization of the [[Predator-Prey Models - Lotka-Volterra vs Rosenzweig-MacArthur|Lotka-Volterra]] system around its coexistence fixed point.

## 1. Coexistence Equilibrium

Setting $\frac{dx}{dt} = \alpha x - \beta xy = 0$ and $\frac{dy}{dt} = \delta xy - \gamma y = 0$ (excluding $(0,0)$):

$$
(x^*, y^*) = \left(\frac{\gamma}{\delta},\; \frac{\alpha}{\beta}\right)
$$

## 2. Jacobian Matrix

$$
J(x, y) = \begin{pmatrix} \alpha - \beta y & -\beta x \\ \delta y & \delta x - \gamma \end{pmatrix}
$$

Evaluating at $(x^*, y^*)$ zeroes both diagonal entries:

$$
J(x^*, y^*) = \begin{pmatrix} 0 & -\frac{\beta\gamma}{\delta} \\ \frac{\delta\alpha}{\beta} & 0 \end{pmatrix} = \begin{pmatrix} 0 & -A \\ B & 0 \end{pmatrix}, \qquad A, B > 0
$$

## 3. Eigenvalues and Center Classification

Solving $\det(J - \lambda I) = \lambda^2 + AB = 0$ gives **purely imaginary eigenvalues**:

$$
\lambda = \pm i\sqrt{AB} = \pm i\sqrt{\alpha\gamma}
$$

- $\mathrm{Re}(\lambda) < 0$: stable spiral (added when logistic $-x^2$ prey self-limitation puts a negative term on the top-left diagonal).
- $\mathrm{Re}(\lambda) > 0$: unstable spiral.
- $\mathrm{Re}(\lambda) = 0$: **Center** (neutrally stable eternal orbits around $(x^*, y^*)$).

## Connections

- [[Predator-Prey Models - Lotka-Volterra vs Rosenzweig-MacArthur]]
- [[Bifurcation Theory and Double-Well Free Energy]]
