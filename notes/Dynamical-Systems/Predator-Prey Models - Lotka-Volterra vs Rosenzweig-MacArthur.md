---
title: "Predator-Prey Models - Lotka-Volterra vs Rosenzweig-MacArthur"
domain: dynamical-systems
source_book: "Mathematical Ecology / Nonlinear Dynamics (Off-Topic)"
section: "2-Species ODEs & Bifurcations"
tags:
  - off-topic
  - dynamical-systems
  - predator-prey
  - lotka-volterra
  - rosenzweig-macarthur
  - chaos
archive_source: "md/Gemini-Predator-Prey Models and Outcomes.md"
status: synthesized
---

# Predator-Prey Models - Lotka-Volterra vs Rosenzweig-MacArthur

*(Off-topic exploration outside elasticity, saved from Gemini chat archive.)*

Nonlinear population ODEs rarely admit closed-form $x(t)$ formulas; instead, they are solved via **conserved invariants**, **Jacobian linear stability**, and **numerical integration** (Runge-Kutta).

## 1. Classical Lotka-Volterra Model

Assumes unlimited prey exponential growth and insatiable linear predation:

$$
\frac{dx}{dt} = \alpha x - \beta x y, \qquad \frac{dy}{dt} = \delta x y - \gamma y
$$

Possesses a conserved quantity $V(x,y) = \delta x - \gamma \ln x + \beta y - \alpha \ln y = C$, producing **neutral closed orbits** whose amplitude depends on initial conditions (see [[Linear Stability Analysis and Jacobian Eigenvalues (Lotka-Volterra)]]).

## 2. Rosenzweig-MacArthur Model

Adds **logistic prey carrying capacity** ($K$) and **Type II saturating functional response** (handling time $h$, conversion efficiency $e$):

$$
\frac{dx}{dt} = r x \left(1 - \frac{x}{K}\right) - \frac{a x y}{1 + h x}, \qquad \frac{dy}{dt} = e\frac{a x y}{1 + h x} - m y
$$

## 3. Dynamical Outcomes ("Resolutions")

| Outcome | Model Required | Mechanism |
| --- | --- | --- |
| **Neutral Cycles** | Lotka-Volterra | Conserved $V(x,y)$; perturbations shift to a new orbit. |
| **Asymptotic Steady State** | Rosenzweig-MacArthur (low $K$) | Damped oscillations settle to $(x^*, y^*)$. |
| **Stable Limit Cycle** | Rosenzweig-MacArthur (high $K$) | **Hopf Bifurcation** (**Paradox of Enrichment**: increasing $K$ destabilizes equilibrium into self-sustaining cycles). |
| **Finite-Time Extinction** | Stochastic / Allee Effect | Heteroclinic orbit swings $< 1$ individual, or Allee reproduction threshold triggers collapse. |
| **Chaos** | 3-Species / Seasonality / Discrete | Impossible in autonomous 2D continuous flows by **Poincaré–Bendixson theorem**; appears in 3-species food chains (Hastings-Powell), seasonal $r(t)$, or discrete maps. |

## Connections

- [[Linear Stability Analysis and Jacobian Eigenvalues (Lotka-Volterra)]]
- [[Bifurcation Theory and Double-Well Free Energy]] (Hopf vs. pitchfork bifurcations)
