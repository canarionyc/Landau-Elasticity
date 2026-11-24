---
title: "Rayleigh Surface Waves and the Cubic Velocity Equation"
domain: elastodynamics
source_book: "Landau & Lifshitz / Lord Rayleigh (1885)"
section: "§24 (Surface Waves)"
tags:
  - elastodynamics
  - rayleigh-waves
  - surface-waves
  - boundary-conditions
archive_source: "md/Gemini-Elastodynamics_ Wave Velocities and Simulation.md"
status: synthesized
---

# Rayleigh Surface Waves and the Cubic Velocity Equation

When P-waves and vertically polarized SV-waves reflect at the Earth's stress-free surface ($z = 0$), they couple into **Rayleigh surface waves**: retrograde elliptical particle motion in the vertical $(x, z)$ plane whose amplitude decays exponentially with depth over roughly one wavelength.

## Derivation Summary

Assume 2D plane-strain potentials decaying in depth ($z > 0$):

$$
\varphi = A e^{-qz} e^{i(kx-\omega t)}, \qquad \psi = B e^{-sz} e^{i(kx-\omega t)}
$$

where $q^2 = k^2(1 - c_R^2/c_p^2)$ and $s^2 = k^2(1 - c_R^2/c_s^2)$. Enforcing stress-free tractions at $z = 0$ ($\sigma_{zz} = 0, \sigma_{xz} = 0$) and setting the $2\times 2$ determinant to zero yields:

$$
\left(2 - \frac{c_R^2}{c_s^2}\right)^2 = 4\sqrt{1 - \frac{c_R^2}{c_p^2}}\sqrt{1 - \frac{c_R^2}{c_s^2}}
$$

Squaring and substituting $\xi = (c_R/c_s)^2$ and $\gamma = c_s/c_p$ gives the **Rayleigh cubic equation**:

$$
\xi^3 - 8\xi^2 + 8(3 - 2\gamma^2)\xi - 16(1 - \gamma^2) = 0, \qquad \gamma^2 = \frac{1-2\nu}{2(1-\nu)}
$$

## Key Result

Because the equation depends only on $\nu$ (no $k$ or $\omega$ dependence), Rayleigh waves in a homogeneous half-space are **non-dispersive** and always slightly slower than shear waves: $c_R \approx 0.87 c_s$ ($\nu = 0$) up to $c_R \approx 0.955 c_s$ ($\nu \to 0.5$).

## Connections

- [[Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)]]
- [[Love Surface Waves and Crustal Waveguide Dispersion]]
- [[Tuned Mass Dampers and Seismic Wave Control]]
