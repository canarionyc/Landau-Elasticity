---
title: "From Love Scale Factors to Metric Tensors and Christoffel Symbols"
domain: tensor-elasticity
source_book: "A. E. H. Love, Treatise on the Mathematical Theory of Elasticity"
section: "Ch. V, Arts. 96-98"
tags:
  - love
  - metric-tensor
  - christoffel-symbols
  - symbolics-jl
  - covariant-derivative
archive_source: "md/Gemini-AEH Love Elasticity Theory_ Classical to Modern.md"
status: synthesized
---

# From Love Scale Factors to Metric Tensors and Christoffel Symbols

In 1906, Love expressed curvilinear line elements via scale factors ($h_1, h_2, h_3$):

$$
ds^2 = \frac{d\alpha^2}{h_1^2} + \frac{d\beta^2}{h_2^2} + \frac{d\gamma^2}{h_3^2}
$$

and re-derived equilibrium equations manually on curvilinear parallelepipeds for each coordinate system.

## Differential Geometry Replacement

Modern tensor calculus replaces manual geometric projections with a **geometry-agnostic engine**:

1. **Metric Tensor ($g_{ij}$)**: Computed from Jacobian $J$ of mapping $q \mapsto x$: $g_{ij} = J^T J$.
2. **Christoffel Symbols of the Second Kind ($\Gamma^k_{ij}$)**:

$$
\Gamma^k_{ij} = \frac{1}{2} g^{kl}\left(\frac{\partial g_{lj}}{\partial x^i} + \frac{\partial g_{li}}{\partial x^j} - \frac{\partial g_{ij}}{\partial x^l}\right)
$$

3. **Kinematics & Equilibrium**:

$$
\varepsilon_{ij} = \frac{1}{2}\left(\nabla_i u_j + \nabla_j u_i\right), \qquad \nabla_i u_j = \frac{\partial u_j}{\partial x^i} - \Gamma^k_{ij} u_k, \qquad \nabla_j \sigma^{ij} + F^i = 0
$$

## Reusable `Symbolics.jl` Christoffel Engine

```julia
using Symbolics, LinearAlgebra

# Function to compute Christoffel Symbols of the Second Kind (Γ^k_ij)
# Γ^k_ij = 0.5 * g^kl * (∂g_lj/∂x_i + ∂g_li/∂x_j - ∂g_ij/∂x_l)
function compute_christoffel(g, g_inv, vars)
    n = length(vars)
    Γ = zeros(Num, n, n, n)
    for k in 1:n, i in 1:n, j in 1:n
        sum_term = 0
        for l in 1:n
            # Partial derivatives of the metric
            dg_lj_i = Symbolics.derivative(g[l, j], vars[i])
            dg_li_j = Symbolics.derivative(g[l, i], vars[j])
            dg_ij_l = Symbolics.derivative(g[i, j], vars[l])
            sum_term += g_inv[k, l] * (dg_lj_i + dg_li_j - dg_ij_l)
        end
        Γ[k, i, j] = 0.5 * sum_term
    end
    return Γ
end
```

Feeding $g_{ij} = \mathrm{diag}(1, r^2, r^2\sin^2\theta)$ into this engine automatically recovers Navier's spherical equilibrium equation $(\lambda+2\mu)(U'' + \frac{2}{r}U' - \frac{2U}{r^2}) = 0$ from [[Spherically Symmetric Elasticity (Landau §7)]].

## Connections

- [[Modernizing Love's Elasticity with 4th-Order Tensors]]
- [[Spherically Symmetric Elasticity (Landau §7)]]
- [[Toroidal Coordinates and the Grad-Shafranov Equation]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
- [[Finite Elasticity and Anisotropic Curved Solids]]
