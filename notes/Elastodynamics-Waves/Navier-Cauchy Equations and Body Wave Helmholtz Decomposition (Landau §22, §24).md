---
title: "Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)"
domain: elastodynamics
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§22, §24"
tags:
  - landau
  - elastodynamics
  - navier-cauchy
  - p-waves
  - s-waves
  - helmholtz-decomposition
archive_source: "md/Gemini-Elastodynamics_ Wave Velocities and Simulation.md"
status: synthesized
---

# Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)

In a homogeneous, linear, isotropic 3D elastic solid, momentum balance yields the **Navier-Cauchy equation**:

$$
\rho \frac{\partial^2 \mathbf{u}}{\partial t^2} = (\lambda + \mu)\nabla(\nabla \cdot \mathbf{u}) + \mu\nabla^2 \mathbf{u} + \mathbf{f}
$$

## Helmholtz Decomposition: P-Waves and S-Waves

Decomposing $\mathbf{u} = \nabla\varphi + \nabla\times\boldsymbol{\psi}$ decouples the homogeneous ($\mathbf{f}=0$) PDE into two non-dispersive body wave equations:

1. **Primary (P) Compressional Waves** (irrotational, longitudinal):

$$
\frac{\partial^2 \varphi}{\partial t^2} = c_p^2 \nabla^2\varphi, \qquad c_p = \sqrt{\frac{\lambda + 2\mu}{\rho}} = \sqrt{\frac{K + \frac{4}{3}\mu}{\rho}}
$$

2. **Secondary (S) Shear Waves** (equivoluminal, transverse):

$$
\frac{\partial^2 \boldsymbol{\psi}}{\partial t^2} = c_s^2 \nabla^2\boldsymbol{\psi}, \qquad c_s = \sqrt{\frac{\mu}{\rho}}
$$

Fluids cannot sustain shear ($\mu = 0$), so S-waves cannot propagate through the liquid outer core.

## Dependence on Poisson's Ratio ($\nu$)

$$
\frac{c_p}{c_s} = \sqrt{\frac{2(1-\nu)}{1-2\nu}}
$$

- **Cork ($\nu \approx 0$)**: $c_p/c_s = \sqrt{2} \approx 1.41$.
- **Steel / Concrete ($\nu \approx 0.3$)**: $c_p/c_s \approx 1.87$.
- **Incompressible Rubber ($\nu \to 0.5$)**: $c_p/c_s \to \infty$.

## Viscosity & Seismic Timescales

At earthquake frequencies ($0.1\text{–}10\text{ Hz}$), the Deborah number is high: mantle rock behaves as an elastic solid, and viscoelasticity contributes attenuation/damping rather than wave propagation speed. Contrast this 3D non-dispersive behavior with 2D bending dispersion in [[Membrane vs Circular Plate Vibrations and Dispersion]].

## Connections

- [[Rayleigh Surface Waves and the Cubic Velocity Equation]]
- [[Love Surface Waves and Crustal Waveguide Dispersion]]
- [[Tuned Mass Dampers and Seismic Wave Control]]
- [[Membrane vs Circular Plate Vibrations and Dispersion]]
- [[Time vs Space Notation Rosetta Stone and Landau Wavevector k]]
