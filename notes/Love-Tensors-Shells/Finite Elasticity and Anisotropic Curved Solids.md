---
title: "Finite Elasticity and Anisotropic Curved Solids"
domain: tensor-elasticity
source_book: "A. E. H. Love / Continuum Mechanics"
section: "Nonlinear & Aeolotropic Elasticity"
tags:
  - finite-elasticity
  - neo-hookean
  - green-lagrange-strain
  - anisotropy
  - biomaterials
archive_source: "md/Gemini-AEH Love Elasticity Theory_ Classical to Modern.md"
status: synthesized
---

# Finite Elasticity and Anisotropic Curved Solids

Where 2D complex potentials ([[Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities]]) fail—large deformations and 3D curved anisotropy—the metric-tensor formalism extends seamlessly without rewriting the physics engine.

## 1. Finite Elasticity (Reference $G_{ij}$ vs. Current $g_{ij}$)

For large stretches (rubber, biological tissues), linear strain breaks down. Strain is defined directly as the difference between the deformed metric $g_{ij}$ and reference metric $G_{ij}$ (**Green-Lagrange Strain**):

$$
E_{ij} = \frac{1}{2}(g_{ij} - G_{ij})
$$

For an incompressible Neo-Hookean rubber sphere ($R \mapsto r(R) = R\lambda(R)$, $\det F = 1 \implies \lambda_r \lambda_\theta^2 = 1$), principal stretches squared are simply ratios of diagonal metric components $\lambda_i^2 = g_{ii}/G_{ii}$, giving first invariant $I_1 = \lambda_r^2 + 2\lambda_\theta^2$ and strain energy $W = C_1(I_1 - 3)$.

## 2. Anisotropic Materials in Curved Geometries

In curved anisotropic structures (wood grain around a trunk, filament-wound carbon fiber pressure vessels), the stiffness tensor $C^{ijkl}(x)$ rotates with the local curvilinear basis. Writing equilibrium in covariant form $\nabla_j \sigma^{ij} = 0$ with $\sigma^{ij} = C^{ijkl}(x)\varepsilon_{kl}$ automatically incorporates the spatial rotation of material fibers via the Christoffel symbols.

## Connections

- [[From Love Scale Factors to Metric Tensors and Christoffel Symbols]]
- [[Modernizing Love's Elasticity with 4th-Order Tensors]]
- [[Spherically Symmetric Elasticity (Landau §7)]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
