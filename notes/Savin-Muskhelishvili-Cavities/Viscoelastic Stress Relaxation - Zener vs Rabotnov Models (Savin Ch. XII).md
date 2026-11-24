---
title: "Viscoelastic Stress Relaxation - Zener vs Rabotnov Models (Savin Ch. XII)"
domain: viscoelasticity
source_book: "G. N. Savin, Stress Concentration Around Holes"
section: "Chapter XII, §1 (Eqs. XII.15, XII.23, XII.24)"
tags:
  - savin
  - viscoelasticity
  - zener-model
  - rabotnov-model
  - mittag-leffler
  - stress-concentration
archive_source: "md/Gemini-Savin's Stress Analysis.md"
status: synthesized
---

# Viscoelastic Stress Relaxation - Zener vs Rabotnov Models (Savin Ch. XII)

Applying the **Volterra Correspondence Principle** (Savin Ch. XII, §1) to a fixed-displacement hole problem (e.g., a tightened bolt hole in a viscoelastic plate/gasket) replaces Young's modulus with a time-dependent relaxation modulus $E(t)$:

$$
K_t(t) = K_{\text{elastic}} \cdot \frac{E(t)}{E_0}
$$

## 1. Zener Model (Standard Linear Solid — Exponential)

$$
E(t) = E_\infty + (E_0 - E_\infty)e^{-t/\tau}
$$

For a circular hole ($K_{\text{elastic}} = 3.0$), $K_t(t)$ decays exponentially from $3.0$ at $t=0$ toward $3.0(E_\infty/E_0)$ as $t \to \infty$, explaining bolt-preload relaxation in gaskets.

## 2. Rabotnov Fractional Exponential Model (Mittag-Leffler)

Real alloys and polymers relax **faster than exponential** at $t \to 0$ (weakly singular kernel $t^\alpha$, $\alpha \in [-0.9, -0.7]$) and **slower than exponential** at long times (**heavy tail / memory**). Replacing $e^{-t/\tau}$ with the Mittag-Leffler function $E_\lambda(z) = \sum_{k=0}^\infty \frac{z^k}{\Gamma(\lambda k + 1)}$ ($\lambda = 1 + \alpha$):

$$
E(t) = E_\infty + (E_0 - E_\infty) E_{1+\alpha}\!\left(-\beta t^{1+\alpha}\right)
$$

## Python Comparison Implementation

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma

class RabotnovSolver:
    def __init__(self, E_inst=200e3, E_long=100e3, beta=0.5, alpha=-0.7):
        """
        beta: The 'rate' parameter (analogous to 1/tau)
        alpha: The fractional parameter (Savin uses values like -0.7)
               Note: The power of t will be lambda = 1 + alpha.
        """
        self.E0 = E_inst
        self.E_inf = E_long
        self.beta = beta
        self.alpha = alpha
        self.lam = 1 + alpha # The actual power, e.g., 0.3

    def mittag_leffler(self, z, n_terms=100):
        """
        Computes the Mittag-Leffler function E_lambda(z) via series expansion.
        E_lam(z) = sum( z^k / Gamma(lam*k + 1) )
        """
        result = np.zeros_like(z, dtype=float)
        for k in range(n_terms):
            term = (z ** k) / gamma(self.lam * k + 1)
            result += term
        return result

    def calculate_kt_history(self, time_array):
        """
        Calculates stress concentration decay using Rabotnov kernel.
        Argument for ML function is -beta * t^(1+alpha)
        """
        z_arg = -self.beta * (time_array ** self.lam)
        decay_factor = self.mittag_leffler(z_arg)
        E_t = self.E_inf + (self.E0 - self.E_inf) * decay_factor
        kt_history = 3.0 * (E_t / self.E0)
        return kt_history
```

## Connections

- [[Elasto-Plastic Slip Zones and Ductile Fracture Limits (Savin Ch. V & VIII)]]
- [[Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole]]
- [[Conformal Mapping and Muskhelishvili Complex Potentials for Elliptical Cavities]]
