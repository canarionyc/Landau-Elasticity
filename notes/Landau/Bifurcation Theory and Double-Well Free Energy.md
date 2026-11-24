---
title: "Bifurcation Theory and Double-Well Free Energy"
domain: physics-and-elasticity
source_book: "Landau & Lifshitz (Phase Transitions & Elasticity)"
section: "General Stability / Symmetry Breaking"
tags:
  - landau
  - bifurcation
  - symmetry-breaking
  - free-energy
  - catastrophe-theory
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Bifurcation Theory and Double-Well Free Energy

The transition of a buckled rod from a straight unstable equilibrium to a bent stable state is an instance of **spontaneous symmetry breaking** (pitchfork bifurcation). Near the critical threshold $P_c$, diverse physical systems obey the universal **Landau Free Energy** expansion:

$$
F(x) = \alpha(P_c - P)x^2 + \beta x^4
$$

- **$P < P_c$**: The $x^2$ coefficient is positive $\implies$ single **U-shaped well**, symmetric state $x = 0$ is stable.
- **$P > P_c$**: The $x^2$ coefficient turns negative $\implies$ **W-shaped double well**; $x = 0$ becomes an unstable local maximum and the system settles into one of two broken-symmetry minima $\pm x_0$.

## Analogous Physical Systems

1. **Ferromagnetism ("Mexican Hat" Potential)**: $F(M) = a(T)M^2 + bM^4$. Below the Curie temperature ($T < T_c$), $a(T) < 0$ and the iron block spontaneously magnetizes along a chosen direction.
2. **Rotating Fluid Planets (Maclaurin $\to$ Jacobi Ellipsoids)**: Above a critical angular momentum $L$ (*Landau & Lifshitz Vol. 6*), the axisymmetric oblate Maclaurin spheroid becomes unstable and breaks rotational symmetry into a triaxial Jacobi ellipsoid.
3. **Rayleigh-Bénard Convection**: Above a critical temperature gradient $\Delta T$, quiescent conductive fluid destabilizes and self-organizes into convective rolls.
4. **Cusp Catastrophe (Snap-Through Buckling)**: Unlike the gradual pitchfork of a column, a shallow arch or dome under downward load loses an energy barrier abruptly and snaps violently to an inverted equilibrium.

## Connections

- [[Euler Elastica and Post-Buckling Stability]]
- [[Predator-Prey Models - Lotka-Volterra vs Rosenzweig-MacArthur]] (Hopf bifurcation in dynamical systems)
- [[Helicoidal Rods, Coil Springs, and Elastic Free Energy (Landau §19)]]
