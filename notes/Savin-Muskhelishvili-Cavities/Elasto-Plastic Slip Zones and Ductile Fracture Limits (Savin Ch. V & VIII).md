---
title: "Elasto-Plastic Slip Zones and Ductile Fracture Limits (Savin Ch. V & VIII)"
domain: fracture-mechanics
source_book: "G. N. Savin, Stress Concentration Around Holes"
section: "Ch. V §8, Ch. VIII §2-3 (Eq. VIII.20)"
tags:
  - savin
  - fracture-mechanics
  - elasto-plasticity
  - dugdale-barenblatt
  - slip-zone
archive_source: "md/Gemini-Savin's Stress Analysis.md"
status: synthesized
---

# Elasto-Plastic Slip Zones and Ductile Fracture Limits (Savin Ch. V & VIII)

In semi-plastic metals such as steel, whether a fracture remains governed by linear elastic stress concentration ("escapes plasticity" / small-scale yielding) or becomes dominated by yielding ("trapped in plasticity" / large-scale yielding) depends on the **Slip Zone** size relative to the crack half-length $l$.

## Dugdale-Barenblatt Slip Zone Model (Eq. VIII.20)

Modeling the plastic zone ahead of a crack of half-length $l$ as a thin collinear slip strip extending to effective half-length $a = l + \epsilon$ under applied tension $p$ and yield strength $\sigma_0$:

$$
p = \frac{2\sigma_0}{\pi}\arccos\left(\frac{l}{a}\right) \iff a = \frac{l}{\cos\left(\frac{\pi p}{2\sigma_0}\right)}
$$

## Regimes & Macroscopic Breakdown Condition

1. **Small-Scale Yielding ($p \ll \sigma_0$, $\epsilon \ll l$)**: $\cos\left(\frac{\pi p}{2\sigma_0}\right) \approx 1 \implies a \approx l$. Linear elastic fracture mechanics ($K_I$, Griffith, Savin Ch. II) remains valid.
2. **Large-Scale / General Yielding ($p \to \sigma_0$)**: $\cos\left(\frac{\pi p}{2\sigma_0}\right) \to 0 \implies a \to \infty$. Once the plastic zone exceeds roughly $10\%$ of crack length ($\epsilon \not\ll l$), linear elasticity fails and elasto-plastic equations (Savin Ch. V) are required.

## Connections

- [[Viscoelastic Stress Relaxation - Zener vs Rabotnov Models (Savin Ch. XII)]]
- [[Hierarchy of Plane-Stress Elliptical Cavity Boundary Conditions]]
- [[Kirsch Circle, Harmonic Ellipse, and Mansfield Neutral Hole]]
