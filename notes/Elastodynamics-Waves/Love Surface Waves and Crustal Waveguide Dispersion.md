---
title: "Love Surface Waves and Crustal Waveguide Dispersion"
domain: elastodynamics
source_book: "A. E. H. Love (1911) / Seismology"
section: "Surface SH Waveguide"
tags:
  - elastodynamics
  - love-waves
  - dispersion
  - seismology
  - specfem3d
archive_source: "md/Gemini-Elastodynamics_ Wave Velocities and Simulation.md"
status: synthesized
---

# Love Surface Waves and Crustal Waveguide Dispersion

Unlike [[Rayleigh Surface Waves and the Cubic Velocity Equation|Rayleigh waves]], **Love waves** (horizontally polarized SH-waves, $v(x,z,t)$ with zero vertical motion) **cannot exist in a homogeneous half-space**. They require a slower crustal waveguide layer (thickness $H$, shear speed $\beta_1$, modulus $\mu_1$) over a faster mantle half-space ($\beta_2 > \beta_1$, modulus $\mu_2$).

## Transcendental Dispersion Relation

Trapping oscillatory SH modes in $0 \le z \le H$ ($\beta_1 < c < \beta_2$) and evanescent decay in $z > H$, then matching:
1. Stress-free top surface ($\partial_z V_1(0) = 0$),
2. Continuous displacement at $z = H$ ($V_1(H) = V_2(H)$),
3. Continuous shear stress at $z = H$ ($\mu_1 V'_1(H) = \mu_2 V'_2(H)$),

gives the **Love wave dispersion relation**:

$$
\tan\left(kH\sqrt{\frac{c^2}{\beta_1^2} - 1}\right) = \frac{\mu_2 \sqrt{1 - \frac{c^2}{\beta_2^2}}}{\mu_1 \sqrt{\frac{c^2}{\beta_1^2} - 1}}
$$

## Physical Consequence & Simulation Standards

- **Dispersion**: Low-frequency (long-wavelength) modes sample the faster mantle ($c \to \beta_2$); high-frequency modes stay trapped in the slower crust ($c \to \beta_1$). Seismic impulses stretch into prolonged horizontal shaking that sweeps across building resonant frequencies (see [[Tuned Mass Dampers and Seismic Wave Control]]).
- **3D Open-Source Solvers**: Pure elastodynamic wavefields (P, S, Rayleigh, Love) are simulated with **SPECFEM3D** (Spectral Element Method) and **SeisSol** (ADER-DG), rather than CFD-focused OpenFOAM.

## Connections

- [[Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)]]
- [[Rayleigh Surface Waves and the Cubic Velocity Equation]]
- [[Tuned Mass Dampers and Seismic Wave Control]]
- [[Membrane vs Circular Plate Vibrations and Dispersion]]
