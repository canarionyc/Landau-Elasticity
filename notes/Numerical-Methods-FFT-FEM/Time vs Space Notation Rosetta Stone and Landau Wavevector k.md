---
title: "Time vs Space Notation Rosetta Stone and Landau Wavevector k"
domain: numerical-methods
source_book: "Landau & Lifshitz (Vols. 2, 3, 7) / Spectral Computing"
section: "Fourier Conventions"
tags:
  - fft
  - landau
  - wavenumber
  - notation
  - conventions
archive_source: "md/Gemini-FFT_ Engineering and Implementation.md"
status: synthesized
---

# Time vs Space Notation Rosetta Stone and Landau Wavevector k

Mixing cycle-based frequencies and radian-based wavenumbers is the #1 source of off-by-$2\pi$ errors in spectral PDE solvers.

## Time Domain vs. Space Domain Rosetta Stone

| Concept | Temporal Domain | Spatial Domain | SI Unit |
| --- | --- | --- | --- |
| **Independent Variable** | $t$ (time) | $x$ (position) | $\text{m}$ |
| **Fundamental Period** | $T$ (period) | $\lambda$ (wavelength) | $\text{m}$ |
| **Cycles Frequency** | $\nu$ or $f$ ($\text{Hz}$) | $\sigma$ or $\xi$ ($\tilde{\nu}$ in spectroscopy) | $\text{cycles/m}$ ($\text{m}^{-1}$) |
| **Angular Frequency** | $\omega = 2\pi f$ | $k = 2\pi \sigma = \frac{2\pi}{\lambda}$ | $\text{rad/m}$ ($\text{m}^{-1}$) |

## Library Gotcha: `fftfreq` Returns Cycles, Not $k$

Both `scipy.fft.fftfreq` and `FFTW.fftfreq` return **cycles per unit length** ($\sigma$). You must multiply by $2\pi$ before spectral differentiation:

```julia
# In Julia
freqs_cycles = rfftfreq(N, 1/N)  # Cycles per unit length
k_landau = 2π * freqs_cycles     # The 'k' used in L&L books
```

## Landau & Lifshitz Convention

Across the *Course of Theoretical Physics* (Vol. 2 *Classical Theory of Fields*, Vol. 3 *Quantum Mechanics*, Vol. 7 *Theory of Elasticity*), Landau exclusively uses **angular wavenumber $k$** and the phase convention:

$$
e^{i(\mathbf{k}\cdot \mathbf{r} - \omega t)}, \qquad p = \hbar k, \qquad k^\mu = \left(\frac{\omega}{c}, \mathbf{k}\right), \qquad k_\mu x^\mu = \omega t - \mathbf{k}\cdot\mathbf{r}
$$

When matching FFT signs to L&L wave equations, verify forward/inverse transform sign conventions so $+i\mathbf{k}\cdot\mathbf{r}$ and $-i\omega t$ do not propagate waves backward in time.

## Connections

- [[Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev]]
- [[Navier-Cauchy Equations and Body Wave Helmholtz Decomposition (Landau §22, §24)]]
- [[Rayleigh Surface Waves and the Cubic Velocity Equation]]
