---
title: "Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev"
domain: numerical-methods
source_book: "Spectral Methods / FFTW.jl"
section: "Spatial PDEs & Boundary Conditions"
tags:
  - fft
  - spectral-methods
  - boundary-conditions
  - julia
  - fftw
  - chebyshev
archive_source: "md/Gemini-FFT_ Engineering and Implementation.md"
status: synthesized
---

# Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev

Using FFTs for spatial PDEs ($\nabla^2 u = f$) replaces finite-difference stencils with algebraic multiplication in wavenumber space ($\frac{d}{dx} \leftrightarrow ik$, $\frac{d^2}{dx^2} \leftrightarrow -k^2$), achieving **exponential (spectral) convergence**.

## The Boundary Condition Matching Table

Standard FFT assumes **periodic** boundaries ($u(0)=u(L), u'(0)=u'(L)$); applying it to non-periodic problems causes **Gibbs ringing**. Instead, choose the transform whose implied symmetry matches the physical BC:

| Transform | Implied Symmetry | Boundary Condition | Julia `FFTW.jl` / Method |
| --- | --- | --- | --- |
| **Standard FFT** | Periodic | $u(0)=u(L),\; u'(0)=u'(L)$ | `fft`, `rfft` |
| **DST-I** | Odd around endpoints | **Dirichlet**: $u(0)=0,\; u(L)=0$ | `FFTW.r2r(..., FFTW.RODFT00)` (interior nodes only; norm $2(N+1)$) |
| **DCT-I** | Even around endpoints | **Neumann**: $u'(0)=0,\; u'(L)=0$ | `FFTW.r2r(..., FFTW.REDFT00)` (remove $k=0$ mean singularity; norm $2(N-1)$) |
| **Chebyshev Matrix** | Non-periodic Gauss-Lobatto | **Robin / Ritz / Mixed**: $u + \alpha u' = 0$ | Row replacement on $L_{op} = -D^2$ |

## Engineering & HPC Best Practices

- **FFT Homogenization**: Periodic RVEs of composite microstructures solve the Lippmann-Schwinger equation iteratively between strain and stress spaces faster than FEM.
- **Dealiasing**: Nonlinear products in spectral space generate ghost high frequencies; zero-pad using the **3/2 rule** (or 2/3 truncation).
- **Memory Layout & Stride**: Transpose 2D/3D arrays (slab/pencil decomposition via `Pencils.jl`) so transformed axes stay contiguous in column-major (Julia/Fortran) or row-major (C/Python) memory.

## Julia Spectral Poisson Solvers (`-u''(x) = f(x)`)

```julia
using FFTW
using LinearAlgebra

# CASE 1: DIRICHLET BCs (u=0 at ends) - Solver: Discrete Sine Transform (DST-I)
function solve_dirichlet(f, x)
    f_vals = f.(x)
    f_hat = FFTW.r2r(f_vals, FFTW.RODFT00)
    k = 1:length(x)
    eigenvalues = k.^2
    u_hat = f_hat ./ eigenvalues
    u_sol = FFTW.r2r(u_hat, FFTW.RODFT00) / (2 * (length(x) + 1))
    return [0; u_sol; 0] # Pad with boundary zeros
end

# CASE 2: NEUMANN BCs (u'=0 at ends) - Solver: Discrete Cosine Transform (DCT-I)
function solve_neumann(f, x)
    f_vals = f.(x)
    f_vals = f_vals .- mean(f_vals)
    f_hat = FFTW.r2r(f_vals, FFTW.REDFT00)
    k = 0:(length(x)-1)
    eigenvalues = k.^2
    eigenvalues[1] = 1.0 # Avoid Inf
    u_hat = f_hat ./ eigenvalues
    u_hat[1] = 0.0       # Enforce 0 mean
    u_sol = FFTW.r2r(u_hat, FFTW.REDFT00) / (2 * (length(x) - 1))
    return u_sol
end
```

*(Robin/Chebyshev collocation replaces boundary rows `L_op[1, :]` and `L_op[end, :]` of $-D^2$ with $u \pm u' = 0$.)*

## Connections

- [[Time vs Space Notation Rosetta Stone and Landau Wavevector k]]
- [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
- [[Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)]] (Robin BC)
