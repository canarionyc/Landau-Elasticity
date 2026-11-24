---
title: "Euler Elastica and Post-Buckling Stability"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§19, §21 (Large Deflections)"
tags:
  - landau
  - elastica
  - nonlinear-elasticity
  - elliptic-integrals
  - buckling
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Euler Elastica and Post-Buckling Stability

Linear buckling assumes small slopes ($w' \ll 1$), leaving the post-buckled amplitude indeterminate at $P = P_{cr}$ and predicting collapse to infinite deflection for $P > P_{cr}$. Retaining exact curvature $\frac{d\theta}{ds}$ yields the nonlinear pendulum-like **Elastica** equation:

$$
EI \frac{d^2\theta}{ds^2} + P \sin\theta = 0
$$

## Exact Elliptic Integral Solution

Parametrizing by the maximum inflection angle $\alpha$ (with elliptic modulus $p = \sin(\alpha/2)$, $m = p^2$), the load required to hold a pinned-pinned column of length $L$ at angle $\alpha$ increases monotonically above $P_{cr}$:

$$
\frac{P}{P_{cr}} = \left(\frac{2 K(m)}{\pi}\right)^2
$$

where $K(m) = \int_0^{\pi/2} \frac{d\theta}{\sqrt{1 - m\sin^2\theta}}$ is the complete elliptic integral of the first kind. Because the vertical projected span shortens ($L_{proj} = L \left(\frac{2 E(m)}{K(m)} - 1\right) < L$), the destabilizing moment arm shrinks, making post-buckling **supercritically stable**.

## Free Energy Landscape (Double Well)

Total potential energy $\Pi = U - W$ compares bending strain energy $U = \frac{1}{2} EI \int_0^L (\theta')^2\,ds$ against load work $W = P(L - L_{proj})$:
- **Linear approx**: $\Pi$ is an inverted parabola for $P > P_{cr}$ (unstable collapse).
- **Exact Elastica**: $\Pi(\alpha)$ bends upward at large $\alpha$, forming a **double-well ("W") potential** with two symmetric stable minima at $\pm \alpha_{eq}$.

## Julia Implementation (Elliptic Integral Solver)

```julia
using Plots
using QuadGK # Standard library for numerical integration

# 1. Helper: Complete Elliptic Integral of the First Kind K(m)
# We compute this numerically to avoid external heavy dependencies
function EllipticK(m)
    # K(m) = integral from 0 to pi/2 of dtheta / sqrt(1 - m*sin(theta)^2)
    f(theta) = 1.0 / sqrt(1.0 - m * sin(theta)^2)
    val, err = quadgk(f, 0.0, pi/2)
    return val
end

function solve_elastica_shape(alpha, n_points=100)
    p = sin(alpha / 2.0)
    m = p^2
    K_val = EllipticK(m)
    load_ratio = (2.0 * K_val / pi)^2

    theta_steps = range(0, alpha - 1e-5, length=n_points ÷ 2)
    x_half = Float64[0.0]
    y_half = Float64[0.0]
    current_x = 0.0
    current_y = 0.0

    for i in 1:length(theta_steps)-1
        th_mid = (theta_steps[i] + theta_steps[i+1])/2
        d_th = theta_steps[i+1] - theta_steps[i]
        denom = sqrt(sin(alpha/2)^2 - sin(th_mid/2)^2)
        ds_dth = 0.5 / denom # The factor 1/2k is absorbed in scale
        current_x += cos(th_mid) * ds_dth * d_th
        current_y += sin(th_mid) * ds_dth * d_th
        push!(x_half, current_x)
        push!(y_half, current_y)
    end

    x_half .*= 1.0 / K_val
    y_half .*= 1.0 / K_val

    X_deflection = [reverse(y_half); -y_half]
    Y_vertical   = [0.5 .- reverse(x_half); 0.5 .+ x_half]
    return X_deflection, Y_vertical, load_ratio
end
```

## Connections

- [[Small Rod Deflection and Euler Buckling (Landau §20-§21)]]
- [[Bifurcation Theory and Double-Well Free Energy]]
- [[Modernizing Love's Elasticity with 4th-Order Tensors]] (Love Ch. XVIII–XIX Elastica)
