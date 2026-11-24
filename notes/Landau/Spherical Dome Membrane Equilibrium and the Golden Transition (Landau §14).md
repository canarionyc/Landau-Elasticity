---
title: "Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)"
domain: shells-and-architecture
source_book: "Landau & Lifshitz, Theory of Elasticity"
section: "§14"
tags:
  - landau
  - shells
  - dome
  - golden-ratio
  - robin-bc
  - rhino-grasshopper
archive_source: "md/Gemini-Landau Elasticity Examples in Julia (2).md"
status: synthesized
---

# Spherical Dome Membrane Equilibrium and the Golden Transition (Landau §14)

Unlike flat plates that carry loads by bending, curved shells carry self-weight $q$ (weight per unit area) primarily via **in-plane membrane forces**: meridional $N_\theta$ and hoop $N_\phi$.

## Membrane Equations for a Spherical Dome

For a sphere of radius $a$ and zenith angle $\theta$ ($\theta = 0$ at crown, $\theta = 90^\circ$ at base):

1. **Meridional Equilibrium** (cap weight $W = 2\pi a^2 q(1-\cos\theta)$ balanced around $2\pi a \sin\theta$):

$$
N_\theta = -\frac{aq}{1 + \cos\theta}
$$

$N_\theta < 0$ everywhere: meridional arches remain in pure compression.

2. **Normal Equilibrium** ($\frac{N_\theta}{a} + \frac{N_\phi}{a} = -q\cos\theta$):

$$
N_\phi = aq\left(\frac{1}{1+\cos\theta} - \cos\theta\right)
$$

## The Golden Transition ($\theta_{cr} \approx 51.8^\circ$)

Setting $N_\phi = 0$ locates where horizontal masonry rings switch from compression to tension:

$$
\cos^2\theta + \cos\theta - 1 = 0 \implies \cos\theta_{cr} = \frac{\sqrt{5}-1}{2} \approx 0.61803 \implies \theta_{cr} \approx 51.8^\circ
$$

Below $51.8^\circ$, unreinforced masonry cracks radially unless restrained by a tension ring—such as Sir Christopher Wren's wrought-iron chain around St. Paul's Cathedral.

## Wren's Iron Chain as a Robin Boundary Condition

A chain of area $A$ and modulus $E_c$ at radius $R$ acts as a radial elastic spring ($N_r - k u_r = 0$ with $k = \frac{E_c A}{R^2}$), interpolating between free ($A=0$) and pinned ($A\to\infty$) edges. This incompatibility between membrane expansion and chain restraint excites the [[Shell Edge Bending Boundary Layer (Geckeler Approximation)]].

## Rhino / Grasshopper 3D Visualization Snippet

```python
#%%
import Rhino.Geometry as rg
import System.Drawing.Color as col
import math

# =============================================================================
# Grasshopper Python Component: Spherical Dome Membrane Stresses
# Inputs: 
#   Radius (float) - Radius of the dome
#   Resolution (int) - Number of mesh faces
# Outputs: 
#   DomeMesh (Mesh) - Colored 3D mesh of the dome
# =============================================================================

def calculate_hoop_force(theta_rad):
    """
    Calculates the normalized hoop force N_phi for a given zenith angle.
    Positive = Tension, Negative = Compression
    """
    cos_t = math.cos(theta_rad)
    return (1.0 / (1.0 + cos_t)) - cos_t

def generate_stress_dome(radius, res):
    mesh = rg.Mesh()
    for i in range(res + 1):
        theta = (math.pi / 2.0) * (i / float(res))
        n_phi = calculate_hoop_force(theta)
        if n_phi > 0:
            intensity = min(int(255 * (n_phi / 0.5)), 255)
            vert_color = col.FromArgb(255, intensity, 0, 0)
        else:
            intensity = min(int(255 * abs(n_phi)), 255)
            vert_color = col.FromArgb(255, 0, 0, intensity)
        for j in range(res):
            phi = 2.0 * math.pi * (j / float(res))
            x = radius * math.sin(theta) * math.cos(phi)
            y = radius * math.sin(theta) * math.sin(phi)
            z = radius * math.cos(theta)
            mesh.Vertices.Add(x, y, z)
            mesh.VertexColors.Add(vert_color)
    for i in range(res):
        for j in range(res):
            next_j = (j + 1) % res
            p1 = i * res + j
            p2 = i * res + next_j
            p3 = (i + 1) * res + next_j
            p4 = (i + 1) * res + j
            mesh.Faces.AddFace(p1, p2, p3, p4)
    mesh.Normals.ComputeNormals()
    return mesh
```

## Connections

- [[Shell Edge Bending Boundary Layer (Geckeler Approximation)]]
- [[General Shell Theory and Paraboloid Membrane Solution (Love Ch. XXIV)]]
- [[Navier Plate Bending and Circular Plate Curvatures (Landau §11-§12)]]
- [[Spatial FFT and Boundary Conditions - DST, DCT, and Chebyshev]] (Robin BCs)
