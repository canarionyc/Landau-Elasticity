---
title: "2D Truss and Frame Data Structures in Rust"
domain: structural-engineering
source_book: "Structural Matrix Analysis (Rust Implementation)"
section: "2D Arena Graph & Local Frame Stiffness"
tags:
  - rust
  - fem
  - truss
  - frame
  - structural-engineering
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# 2D Truss and Frame Data Structures in Rust

To avoid borrow-checker friction in Rust graph structures, structural meshes use an **index-based arena pattern**: `Vec<Node>`, `Vec<Element>`, `Vec<Constraint>`, and `Vec<Load>` storing `usize` node indices rather than object references.

## Truss vs. Welded Frame (Table Frame Example)

- **2D Truss**: Pin-jointed axial bars ($2\text{ DoF/node}: u_x, u_y$; element stiffness uses $EA/L$). A rectangular steel table modeled only with pinned truss joints without diagonal bracing is a **mechanism** that collapses sideways.
- **2D / 3D Frame**: Rigid welded joints transmit moments ($3\text{ DoF/node}$ in 2D: $u_x, u_y, \theta_z$; $6\times 6$ local stiffness matrix combining axial $EA/L$ and Euler-Bernoulli cubic bending $12EI/L^3, 6EI/L^2, 4EI/L, 2EI/L$).

## Rust Arena Structs

```rust
// =========================================================
// 2D TRUSS STRUCTURE IN RUST
// Pattern: Index-Based Graph (Arena style)
// =========================================================

/// Represents a geometric point in 2D space
#[derive(Debug, Clone, Copy)]
pub struct Node {
    pub x: f64,
    pub y: f64,
}

/// Represents a physical connection between two nodes
#[derive(Debug, Clone)]
pub struct Element {
    pub node_ids: [usize; 2], // Indices into the Truss.nodes vector
    pub area: f64,            // Cross-sectional Area
    pub modulus: f64,         // Young's Modulus (E)
}

/// Represents a boundary condition (Support)
/// Using booleans: true = fixed (cannot move), false = free
#[derive(Debug, Clone, Copy)]
pub struct Constraint {
    pub node_id: usize,
    pub fix_x: bool,
    pub fix_y: bool,
}

/// Represents an external force applied to a node
#[derive(Debug, Clone, Copy)]
pub struct Load {
    pub node_id: usize,
    pub fx: f64,
    pub fy: f64,
}
```

## Connections

- [[3D Welded Frame Stiffness, I-Beam Orientation, and Global Assembly]]
- [[FDM vs FEM and Hermite Cubic Shape Functions for Buckling]]
- [[Small Rod Deflection and Euler Buckling (Landau §20-§21)]]
