---
title: "Tuned Mass Dampers and Seismic Wave Control"
domain: structural-dynamics
source_book: "Structural Dynamics & Seismic Design"
section: "TMD Phase & Base Isolation"
tags:
  - structural-dynamics
  - tuned-mass-damper
  - love-waves
  - rayleigh-waves
  - base-isolation
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Tuned Mass Dampers and Seismic Wave Control

Buildings possess very high vertical **axial stiffness** in columns but much lower lateral **bending/shear stiffness**. Consequently, horizontal [[Love Surface Waves and Crustal Waveguide Dispersion|Love waves]] and vertical [[Rayleigh Surface Waves and the Cubic Velocity Equation|Rayleigh waves]] require fundamentally different mitigation strategies.

## Horizontal Love Waves $\to$ Tuned Mass Damper (TMD)

A TMD (typically $1\%\text{–}5\%$ of building mass $m_1$) is tuned to the building's fundamental lateral natural frequency ($\omega_2 = \omega_1 = \sqrt{k_1/m_1}$).

- **Why $90^\circ$ Phase Shift is Ideal**: For a driven damped oscillator, phase lag is:

$$
\phi = \arctan\left(\frac{2\zeta (\omega/\omega_n)}{1 - (\omega/\omega_n)^2}\right)
$$

At exact resonance ($\omega/\omega_n = 1$), $\phi = -90^\circ$. A $90^\circ$ lag means the damper's velocity—and thus the restoring/damping force exerted back onto the building—is $180^\circ$ out of phase with the building's lateral velocity, maximizing energy extraction (whereas a $180^\circ$ phase lag would be purely reactive rather than dissipative).

## Vertical Rayleigh Waves $\to$ Base Isolation

TMDs are ineffective against vertical Rayleigh wave uplift and rocking because vertical structural frequencies are much higher and axial column compression cannot be countered by a small roof pendulum. Vertical shock mitigation relies on **base isolation** (elastomeric/lead-rubber bearings) and foundation decoupling.

## Connections

- [[Love Surface Waves and Crustal Waveguide Dispersion]]
- [[Rayleigh Surface Waves and the Cubic Velocity Equation]]
- [[Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)]]
