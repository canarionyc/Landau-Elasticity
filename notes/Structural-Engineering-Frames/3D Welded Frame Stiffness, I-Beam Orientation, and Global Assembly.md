---
title: "3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly"
domain: structural-engineering
source_book: "3D Matrix Structural Analysis (Rust Implementation)"
section: "12-DoF Frame, I-Beam Sections & Rotation"
tags:
  - rust
  - fem
  - 3d-frame
  - i-beam
  - torsion
  - rotation-matrix
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# 3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly

Upgrading [[2D Truss and Frame Data Structures in Rust|2D frames]] to a **3D welded steel table frame** separates `Section` profiles from `FrameElement` geometry (`6 DoF/node`: $u_x, u_y, u_z, \theta_x, \theta_y, \theta_z$; $12\times 12$ element matrices) and adds a beam `roll_angle` $\gamma$.

## I-Beam Local Axes & Section Properties

- **Local $x$**: Longitudinal axis along beam.
- **Local $z$ (Strong Axis)**: Horizontal centroidal axis; flanges take maximum compression/tension ($I_z$ dominates via $h^3$).
- **Local $y$ (Weak Axis)**: Vertical web axis ($I_y$ scales with flange width $b^3$).
- **Torsional Constant ($J$)**: Open I-sections warp under torsion (see [[Torsion of Non-Circular Rods and Prandtl Soap-Film Analogy (Landau §16)]]), giving $J \approx \frac{1}{3}[2 b t_f^3 + (h - t_f)t_w^3] \ll I_y, I_z$.

$$
A = 2 b t_f + (h - 2 t_f)t_w, \qquad I_z = \frac{1}{12}\left[b h^3 - (b - t_w)(h - 2 t_f)^3\right], \qquad I_y = \frac{1}{6}t_f b^3 + \frac{1}{12}(h - 2 t_f)t_w^3
$$

## 3D Rotation Matrix ($T$) & Global Assembly

1. **Local Axis Basis ($3\times 3$ block $R$)**: Let $\mathbf{u} = \Delta\mathbf{x}/L$ (local $x$). Choose reference vector $\mathbf{r} = (0,0,1)$—**except** for vertical columns ($|u_z| \approx 1$), where $\mathbf{r} = (1,0,0)$ avoids a zero cross product. Compute $\mathbf{w} = \widehat{\mathbf{u}\times\mathbf{r}}$, $\mathbf{v} = \mathbf{w}\times\mathbf{u}$, and apply `roll_angle` $\gamma$ around local $x$.
2. **Block-Diagonal Transformation ($12\times 12$)**: $T = \mathrm{diag}(R, R, R, R)$ rotates nodal translations and rotations identically.
3. **Global Element Stiffness & Assembly**:

$$
K_{\text{elem, global}} = T^T K_{\text{local}} T
$$

Sum into global $K_{global}$ using the 12-entry `dof_map` for `[6*n1 .. 6*n1+5, 6*n2 .. 6*n2+5]`.

## Connections

- [[2D Truss and Frame Data Structures in Rust]]
- [[Torsion of Non-Circular Rods and Prandtl Soap-Film Analogy (Landau §16)]]
- [[Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)]] (slab carried by the table frame)
- [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]
