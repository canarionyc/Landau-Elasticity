---
title: "Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)"
domain: elasticity
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§11, §12"
tags:
  - landau
  - plates
  - biharmonic
  - navier-solution
  - bending-moments
  - inflection-point
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)

Thin isotropic plates under transverse load $P(x, y)$ obey the fourth-order **biharmonic equation** (*Landau §11*):

$$
D \Delta^2 \zeta = P(x, y), \qquad D = \frac{Eh^3}{12(1-\sigma^2)}
$$

## Rectangular Simply Supported Slab: Navier Solution (§12)

For a rectangular plate $a \times b$ with $\zeta = 0$ and $M_n = 0$ on all edges under uniform load $q_0$:

$$
\zeta(x, y) = \frac{16 q_0}{\pi^6 D} \sum_{m=1,3,5\dots}^\infty \sum_{n=1,3,5\dots}^\infty \frac{\sin\left(\frac{m\pi x}{a}\right)\sin\left(\frac{n\pi y}{b}\right)}{mn\left[\left(\frac{m}{a}\right)^2 + \left(\frac{n}{b}\right)^2\right]^2}
$$

When $a \gg b$, the $1/b$ terms dominate and the slab carries load almost entirely across the short span as a 1D strip.

## Circular Plate: Moments, Curvatures & Inflection Points (§12)

Even under purely axisymmetric loading, radial bending couples with tangential curvature ($\kappa_\phi = -\frac{1}{r}\frac{d\zeta}{dr}$) via Poisson's ratio $\sigma$:

$$
M_r = -D\left(\frac{d^2\zeta}{dr^2} + \frac{\sigma}{r}\frac{d\zeta}{dr}\right), \qquad M_\phi = -D\left(\frac{1}{r}\frac{d\zeta}{dr} + \sigma\frac{d^2\zeta}{dr^2}\right)
$$

### Case A: Clamped Circular Plate ($\zeta(R) = 0, \zeta'(R) = 0$)

$$
\zeta(r) = \frac{P}{64D}(R^2 - r^2)^2
$$

- **Bending Moments**: Positive (bottom tension) at $r = 0$ ($M_r = M_\phi = \frac{PR^2}{16}(1+\sigma)$), reversing to large negative radial moment at the clamped rim ($M_r(R) = -\frac{PR^2}{8}$), explaining circumferential top-surface cracking around column capitals.
- **Curvatures & Inflection Point**: $\kappa_\phi$ remains positive everywhere, whereas $\kappa_r = \frac{P}{16D}(R^2 - 3r^2)$ changes sign at:

$$
r_{\text{inflect}} = \frac{R}{\sqrt{3}} \approx 0.577 R
$$

### Case B: Simply Supported Circular Plate ($\zeta(R) = 0, M_r(R) = 0$)

$$
\zeta(r) = \frac{P(R^2 - r^2)}{64D}\left(\frac{5+\sigma}{1+\sigma}R^2 - r^2\right)
$$

Because $M_r \propto \kappa_r + \sigma \kappa_\phi$ must vanish at $r = R$ while $\kappa_\phi > 0$, $\kappa_r$ still turns slightly negative near the edge! With concrete $\sigma = 0.2$, the inflection point shifts outward to:

$$
r_{\text{inflect}} = R\sqrt{\frac{3+\sigma}{3(1+\sigma)}} \approx 0.942 R
$$

## Connections

- [[Membrane vs Circular Plate Vibrations and Dispersion]]
- [[Plate Buckling and Garland Curves in Box Girder Bridges]]
- [[Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)]]
- [[3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly]] (concrete slab resting on steel table frame)
