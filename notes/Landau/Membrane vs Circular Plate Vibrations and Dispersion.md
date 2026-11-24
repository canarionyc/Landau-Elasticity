---
title: "Membrane vs Circular Plate Vibrations and Dispersion"
domain: elastodynamics-and-plates
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§12, §25 (Plate & Membrane Vibrations)"
tags:
  - landau
  - plates
  - membranes
  - bessel-functions
  - dispersion
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Membrane vs Circular Plate Vibrations and Dispersion

Although visually similar, a tensioned drum membrane and a rigid circular plate differ fundamentally in restoring mechanism, PDE order, boundary conditions, and wave dispersion.

| Property | Drum Membrane | Rigid Plate |
| --- | --- | --- |
| **Restoring Force** | In-plane tension $T$ | Flexural rigidity $D = \frac{Eh^3}{12(1-\sigma^2)}$ |
| **Governing PDE** | $T\nabla^2 z = \rho_{area}\frac{\partial^2 z}{\partial t^2}$ (2nd order) | $-D\nabla^4 z = \rho h \frac{\partial^2 z}{\partial t^2}$ (4th order) |
| **Edge BCs** | 1 BC ($z = 0$) | 2 BCs ($z=0$ and $\partial z/\partial r = 0$ or $M_r = 0$) |
| **Radial Mode Shape** | $J_n(kr)$ | $A J_n(kr) + B I_n(kr)$ |
| **Dispersion** | Non-dispersive ($v_p = v_g = \sqrt{T/\rho_{area}}$) | Dispersive ($v_p = \sqrt{\omega}\left(\frac{D}{\rho h}\right)^{1/4}$, $v_g = 2v_p$) |

## Dispersion Physics

For a bending wave $z = A e^{i(kx - \omega t)}$ in a plate, $D k^4 = \rho h \omega^2$, giving $\omega \propto k^2$:
- **Phase velocity**: $v_p = \omega/k = k\sqrt{D/\rho h} = \sqrt{\omega}(D/\rho h)^{1/4}$.
- **Group velocity**: $v_g = \partial\omega/\partial k = 2 v_p$.

High frequencies outrun low frequencies (metallic "clang" vs. clean drum "thump").

## Characteristic Roots: `scipy.special` vs. Clamped Plate Roots

- **Drum membrane ($J_n(kR) = 0$)**: Tabulated in `scipy.special.jn_zeros(n, k)` (e.g., $J_0$ roots: $2.4048, 5.5201, 8.6537$).
- **Clamped circular plate**: Combining $z(R)=0$ and $z'(R)=0$ for $A J_n(kR) + B I_n(kR)$ yields the coupled characteristic equation:

$$
J_n(\lambda) I'_n(\lambda) - I_n(\lambda) J'_n(\lambda) = 0 \qquad (\lambda = kR)
$$

Because this couples $J_n$ and $I_n$, roots are not in `scipy.special` and must be bracketed numerically (e.g., for $n = 0$: $\lambda_1 \approx 3.196$, $\lambda_2 \approx 6.306$, $\lambda_3 \approx 9.439$).

## Connections

- [[Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)]]
- [[Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)]]
- [[Love Surface Waves and Crustal Waveguide Dispersion]]
