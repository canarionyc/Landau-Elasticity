# Elastodynamics: Wave Velocities and Simulation

**Exported:** 2/22/2026 4:07:28  
**Link:** [https://gemini.google.com/app/7a91529aeeac22a9](https://gemini.google.com/app/7a91529aeeac22a9)

## Prompt:
3D Navier-Cauchy equations of elastodynamics (covered beautifully in Landau & Lifshitz Theory of Elasticity §22 and §24). Hello I was looking into the description of earthquakes and specially on reviewing and updating the landmark collection of books of Landau and Lifshitz. This is a continuation another thread that has grown too large.

## Response:
Gemini said
-----------

The Navier-Cauchy equations serve as the definitive bridge between pure continuum mechanics and the practical realities of structural response. Modernizing the classic, purely analytical approach of Landau and Lifshitz requires contextualizing these equations not just as abstract displacement fields, but as the governing dynamics for seismic energy transfer into the built environment.

Here is a look at the core theory, its decomposition into earthquake body waves, and a computational approach to visualize the material dependencies.

### The Elastodynamic Foundation

For a linear, isotropic, homogeneous elastic medium, the equation of motion derived from the balance of linear momentum and Hooke's Law is the Navier-Cauchy equation:

$$
\rho \frac{\partial ^{2}u}{\partial t^{2}}=\left(\lambda +\mu \right)\nabla \left(\nabla \cdot u\right)+\mu \nabla ^{2}u+f
$$

Where:

*    $u\left(x,t\right)$  is the displacement vector field.
*    $\rho$  is the mass density.
*    $\lambda$  and  $\mu$  are the Lamé parameters ( $\mu$  being the shear modulus).
*    $f$  represents body forces (often neglected when studying wave propagation away from the seismic source).

### Helmholtz Decomposition: The Genesis of Earthquakes

Landau and Lifshitz elegantly demonstrate that an arbitrary vector field can be decomposed into an irrotational (curl-free) scalar potential  $\varphi$  and a solenoidal (divergence-free) vector potential  $\psi$ :

$$
u=\nabla \varphi +\nabla \times \psi
$$

Substituting this into the homogeneous Navier-Cauchy equation ( $f=0$ ) decouples the complex system into two distinct, classic wave equations. This decoupling is the physical basis for the two primary types of earthquake body waves.

**1\. The Primary (P) Wave** Taking the divergence of the Navier-Cauchy equation isolates the volumetric expansion/contraction:

$$
\frac{\partial ^{2}\varphi }{\partial t^{2}}=c_{p2}\nabla ^{2}\varphi wherec_{p}=\sqrt{\frac{\lambda +2\mu }{\rho }}
$$

These are the longitudinal, compressional waves. They arrive first at structural foundations and push/pull the medium in the direction of propagation.

**2\. The Secondary (S) Wave** Taking the curl isolates the equivoluminal shear deformations:

$$
\frac{\partial ^{2}\psi }{\partial t^{2}}=c_{s2}\nabla ^{2}\psi wherec_{s}=\sqrt{\frac{\mu }{\rho }}
$$

These transverse waves arrive second but are responsible for the most severe horizontal accelerations that compromise architectural integrity. Because fluids cannot support shear stress ( $\mu =0$ ), S-waves cannot propagate through the Earth's outer core or local liquid bodies.

* * *

### Computational Update: Wave Velocities and Material Properties

To modernize L&L §22 and §24, we can analyze how the relative speeds of these waves depend entirely on the material's Poisson's ratio ( $\nu$ ). The ratio of the P-wave to S-wave velocity is given by:

$$
\frac{c_{p}}{c_{s}}=\sqrt{\frac{2\left(1-\nu \right)}{1-2\nu }}
$$

Below is a Python snippet to visualize this relationship, structured for cell-based execution and adhering to bilingual rendering requirements.

```
#%%
import numpy as np
import matplotlib.pyplot as plt

# Toggle language: 'EN' or 'ES'
LANG = 'EN'

#%%
# Define Poisson's ratio range (0 to 0.49 to avoid singularity at 0.5)
nu = np.linspace(0, 0.49, 100)
# Calculate velocity ratio cp / cs
vp_vs_ratio = np.sqrt((2 * (1 - nu)) / (1 - 2 * nu))

# Bilingual dictionaries
text = {
    'EN': {
        'title': "Ratio of P-Wave to S-Wave Velocity vs. Poisson's Ratio",
        'xlabel': "Poisson's Ratio ($\\nu$)",
        'ylabel': "Velocity Ratio ($c_p / c_s$)",
        'annot_cork': "Cork ($\\nu \\approx 0$)",
        'annot_steel': "Steel ($\\nu \\approx 0.3$)",
        'annot_rubber': "Rubber ($\\nu \\to 0.5$)"
    },
    'ES': {
        'title': "Razón de Velocidad de Onda P a Onda S vs. Coeficiente de Poisson",
        'xlabel': "Coeficiente de Poisson ($\\nu$)",
        'ylabel': "Razón de Velocidad ($c_p / c_s$)",
        'annot_cork': "Corcho ($\\nu \\approx 0$)",
        'annot_steel': "Acero ($\\nu \\approx 0.3$)",
        'annot_rubber': "Goma ($\\nu \\to 0.5$)"
    }
}

#%%
# Plotting
plt.figure(figsize=(9, 6))
plt.plot(nu, vp_vs_ratio, color='darkred', linewidth=2)

# Annotations for common architectural/structural materials
plt.axvline(x=0.0, color='gray', linestyle='--', alpha=0.5)
plt.text(0.01, 3.5, text[LANG]['annot_cork'], rotation=90)

plt.axvline(x=0.3, color='gray', linestyle='--', alpha=0.5)
plt.text(0.31, 3.5, text[LANG]['annot_steel'], rotation=90)

plt.axvline(x=0.49, color='gray', linestyle='--', alpha=0.5)
plt.text(0.46, 3.5, text[LANG]['annot_rubber'], rotation=90)

plt.title(text[LANG]['title'], fontsize=14, fontweight='bold')
plt.xlabel(text[LANG]['xlabel'], fontsize=12)
plt.ylabel(text[LANG]['ylabel'], fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.ylim(1, 5)
plt.xlim(-0.02, 0.52)
plt.tight_layout()
plt.show()
#%%
```

Would you like to proceed by looking at how to discretize the Navier-Cauchy equations into a 2D Finite Difference Time Domain (FDTD) grid to actually simulate and plot these wavefronts colliding with a solid boundary?

## Prompt:
These are volumetric waves, but when the reach the Earth's crust surface they become the surface waves of Rayleigh and ???

## Response:
Gemini said
-----------

The missing piece of the puzzle is the **Love wave**.

When the volumetric P and S waves (body waves) encounter the Earth's stress-free surface, they reflect, refract, and constructively interfere to form surface waves. Mathematically, this happens because the Navier-Cauchy equations must now satisfy a new boundary condition at the surface ( $z=0$ ): the traction vector must vanish, meaning the stress tensor components  $\sigma _{xz}$ ,  $\sigma _{yz}$ , and  $\sigma _{zz}$  equal zero.

This boundary condition traps seismic energy along the surface, resulting in two distinct types of surface waves that cause the majority of damage to the built environment.

### 1\. Rayleigh Waves (The Ground Roll)

Derived by Lord Rayleigh in 1885, these waves are formed by the coupling and interference of P-waves and SV-waves (vertically polarized shear waves) at the free surface.

*   **Particle Motion:** They induce a retrograde, elliptical particle motion in the vertical plane (parallel to the direction of propagation). As the wave moves forward, the ground is pushed up, backward, down, and forward.
*   **Decay:** The amplitude decays exponentially with depth, typically penetrating only about one wavelength into the medium.
*   **Structural Impact:** For foundations, Rayleigh waves are incredibly destructive because they introduce both violent vertical uplifting and a "pitching" or rocking motion to the structure, forcing columns and load-bearing walls out of plumb.

### 2\. Love Waves (The Horizontal Shear)

Named after A.E.H. Love who mathematically predicted them in 1911, these are purely transverse waves. Interestingly, Love waves cannot exist in a perfectly homogeneous half-space. They _require_ a layered medium—specifically, a superficial layer with a lower shear velocity over a higher-velocity half-space (exactly like the Earth's crust resting on the denser mantle).

*   **Particle Motion:** They are formed by the total internal reflection and constructive interference of SH-waves (horizontally polarized shear waves) trapped within the upper waveguide layer. The particle motion is purely horizontal and perpendicular to the direction of propagation. There is absolutely no vertical displacement.
*   **Velocity:** Love waves generally travel slightly faster than Rayleigh waves but slower than body S-waves ( $c_{R}<c_{L}<c_{S}$ ). They are also highly dispersive, meaning their velocity depends on their frequency.
*   **Structural Impact:** From a structural perspective, Love waves are arguably the most dangerous. They apply massive lateral shearing forces to the base of a structure. Because buildings are typically designed to withstand tremendous vertical static loads (gravity) but are often more vulnerable to severe lateral dynamic loads, the alternating horizontal tearing of a Love wave easily compromises stiff, unreinforced elements.

* * *

The beauty of these waves in elastodynamics is how they emerge purely from the algebraic constraints of boundary conditions applied to the Helmholtz decomposition we looked at earlier.

Would you like to walk through the mathematical derivation of the Rayleigh wave velocity—the classic process of applying the stress-free boundary conditions to yield the famous cubic equation for  $\left(c_{R}/c_{S}\right)^{2}$ ?

## Prompt:
Yes please

## Response:
Gemini said
-----------

This derivation is one of the most elegant intersections of mathematical physics and structural engineering. It takes the abstract field equations of a continuum and, simply by enforcing a boundary condition, predicts the exact ground-rolling phenomena that architects must design foundations to withstand.

Let's walk through the classic derivation of the Rayleigh wave velocity, moving from the wave potentials to the famous cubic equation.

### 1\. The Setup: Potentials and the Plane Wave Ansatz

We consider a semi-infinite elastic half-space (the Earth). Let the free surface be the  $xy$ \-plane ( $z=0$ ), with the  $z$ \-axis pointing positively downwards into the medium. We assume a plane wave propagating in the  $x$ \-direction, making this a 2D plane strain problem where displacements are  $u=\left(u_{x},0,u_{z}\right)$ .

Using the Helmholtz decomposition from earlier, we express the displacements using a scalar potential  $\varphi$  (for P-waves) and a vector potential  $\psi$  (for S-waves). In 2D, the vector potential only needs its  $y$ \-component, which we will simply call  $\psi$ :

$$
u_{x}=\frac{\partial \varphi }{\partial x}-\frac{\partial \psi }{\partial z}
$$

$$
u_{z}=\frac{\partial \varphi }{\partial z}+\frac{\partial \psi }{\partial x}
$$

Because Rayleigh waves are _surface_ waves, their amplitude must decay as depth  $z$  increases. We postulate harmonic solutions that travel in  $x$  but decay exponentially in  $z$ :

$$
\varphi =Ae^{-qz}e^{i\left(kx-\omega t\right)}
$$

$$
\psi =Be^{-sz}e^{i\left(kx-\omega t\right)}
$$

Where:

*    $k$  is the wavenumber and  $\omega$  is the angular frequency. The Rayleigh wave speed is  $c_{R}=\omega /k$ .
*    $q$  and  $s$  are the attenuation constants for the P-wave and S-wave components, respectively.

Substituting these into the decoupled wave equations ( $\nabla ^{2}\varphi =\frac{1}{c_{p2}}\frac{\partial ^{2}\varphi }{\partial t^{2}}$  and  $\nabla ^{2}\psi =\frac{1}{c_{s2}}\frac{\partial ^{2}\psi }{\partial t^{2}}$ ) forces the attenuation constants to be:

$$
q^{2}=k^{2}\left(1-\frac{c_{R2}}{c_{p2}}\right)
$$

$$
s^{2}=k^{2}\left(1-\frac{c_{R2}}{c_{s2}}\right)
$$

### 2\. The Stress-Free Boundary Conditions

At the Earth's surface ( $z=0$ ), there is no material above to push or pull back. Therefore, the normal stress  $\sigma _{zz}$  and the shear stress  $\sigma _{xz}$  must vanish:

$$
\sigma _{zz}​_{z=0}=\lambda \left(\frac{\partial u_{x}}{\partial x}+\frac{\partial u_{z}}{\partial z}\right)+2\mu \frac{\partial u_{z}}{\partial z}=0
$$

$$
\sigma _{xz}​_{z=0}=\mu \left(\frac{\partial u_{x}}{\partial z}+\frac{\partial u_{z}}{\partial x}\right)=0
$$

By substituting our potential-derived displacements ( $u_{x},u_{z}$ ) into these stress equations, and evaluating at  $z=0$ , we obtain a homogenous system of two linear equations for the amplitudes  $A$  and  $B$ :

$$
A\left[2\mu q^{2}+\lambda \left(q^{2}-k^{2}\right)\right]-B\left[2i\mu ks\right]=0
$$

$$
A\left[-2ikq\right]+B\left[s^{2}+k^{2}\right]=0
$$

### 3\. The Rayleigh Equation

For a non-trivial solution to exist (meaning an earthquake actually happens,  $A\ne 0,B\ne 0$ ), the determinant of the coefficient matrix must be zero. After substituting the relationships for  $q^{2}$  and  $s^{2}$ , factoring out  $k$ , and utilizing the identity  $c_{s2}=\mu /\rho$ , the determinant reduces to the foundational Rayleigh equation:

$$
\left(2-\frac{c_{R2}}{c_{s2}}\right)^{2}=4\sqrt{1-\frac{c_{R2}}{c_{p2}}}\sqrt{1-\frac{c_{R2}}{c_{s2}}}
$$

To solve this, we square both sides to remove the radicals. Let  $\xi =\left(c_{R}/c_{s}\right)^{2}$  and let the material parameter be  $\gamma =c_{s}/c_{p}$ . Expanding the terms yields the famous Rayleigh cubic equation in  $\xi$ :

$$
\xi ^{3}-8\xi ^{2}+8\left(3-2\gamma ^{2}\right)\xi -16\left(1-\gamma ^{2}\right)=0
$$

This single algebraic equation entirely governs the speed of the destructive surface roll relative to the shear wave, depending strictly on the material's Poisson's ratio (since  $\gamma ^{2}=\frac{1-2\nu }{2\left(1-\nu \right)}$ ).

### Computational Implementation

Since there isn't a simple, universally neat analytical form for the roots of this cubic across all valid Poisson's ratios, we solve it numerically. Here is the Python implementation to map the relationship, keeping our cell structure and bilingual parameters intact.

```
#%%
import numpy as np
import matplotlib.pyplot as plt

# Toggle language: 'EN' or 'ES'
LANG = 'EN'

#%%
# Define Poisson's ratio range
nu_vals = np.linspace(0, 0.49, 100)
cr_cs_ratio = []

for nu in nu_vals:
    # Material parameter gamma^2 = (cs/cp)^2
    gamma_sq = (1 - 2 * nu) / (2 * (1 - nu))
    
    # Coefficients of the Rayleigh cubic: xi^3 - 8*xi^2 + 8*(3 - 2*gamma^2)*xi - 16*(1 - gamma^2) = 0
    coeffs = [1, -8, 8 * (3 - 2 * gamma_sq), -16 * (1 - gamma_sq)]
    
    # Calculate roots
    roots = np.roots(coeffs)
    
    # The physical root is the real root between 0 and 1 (since c_R < c_S)
    real_roots = np.real(roots[np.iscomplex(roots) == False])
    valid_xi = real_roots[(real_roots > 0) & (real_roots < 1)][0]
    
    # c_R / c_S is the square root of xi
    cr_cs_ratio.append(np.sqrt(valid_xi))

#%%
# Bilingual dictionaries
text = {
    'EN': {
        'title': "Rayleigh Wave Velocity relative to Shear Wave ($c_R/c_S$)",
        'xlabel': "Poisson's Ratio ($\\nu$)",
        'ylabel': "Velocity Ratio ($c_R / c_S$)"
    },
    'ES': {
        'title': "Velocidad de Onda de Rayleigh relativa a la Onda de Corte ($c_R/c_S$)",
        'xlabel': "Coeficiente de Poisson ($\\nu$)",
        'ylabel': "Razón de Velocidad ($c_R / c_S$)"
    }
}

#%%
# Plotting
plt.figure(figsize=(9, 6))
plt.plot(nu_vals, cr_cs_ratio, color='navy', linewidth=2.5)

plt.title(text[LANG]['title'], fontsize=14, fontweight='bold')
plt.xlabel(text[LANG]['xlabel'], fontsize=12)
plt.ylabel(text[LANG]['ylabel'], fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)

# Highlight that Rayleigh waves are always slower than Shear waves
plt.axhline(y=1.0, color='red', linestyle='--', alpha=0.6, label='$c_R = c_S$')
plt.ylim(0.85, 1.05)
plt.xlim(-0.02, 0.52)

if LANG == 'EN':
    plt.legend(loc='lower right')
else:
    plt.legend(['Razón $c_R / c_S$', 'Límite Superior ($c_R = c_S$)'], loc='lower right')

plt.tight_layout()
plt.show()
#%%
```

The resulting plot will explicitly show that  $c_{R}$  is always slightly less than  $c_{S}$  (ranging from about  $0.87c_{s}$  for materials with  $\nu \approx 0$  up to  $0.95c_{s}$  for incompressible materials like rubber).

Would you like to shift from the continuous analytical domain and explore how we would discretize these specific stress-free boundary conditions into a 2D Finite Difference grid to simulate a seismic pulse hitting a building foundation?

## Prompt:
But before there an analytical solution for the Love wave?

## Response:
Gemini said
-----------

Yes, absolutely. The analytical derivation of the Love wave is a perfect follow-up. While the Rayleigh wave emerges from the boundary conditions of a homogeneous half-space, the Love wave mathematically proves why surface geology—specifically, the layering of the Earth's crust—is so critical to earthquake engineering.

To find the analytical solution, we must abandon the homogeneous half-space. Love waves **cannot exist** mathematically without a waveguide: a slower, less dense superficial layer resting on top of a faster, denser half-space.

Here is the derivation of the Love wave dispersion relation, which reveals why these waves are so structurally destructive.

### 1\. The Setup: A Layer Over a Half-Space

Imagine a superficial crustal layer (Layer 1) of thickness  $H$  resting on a semi-infinite mantle (Layer 2). The free surface is at  $z=0$ , and the interface is at  $z=H$ .

*   **Layer 1:** Shear modulus  $\mu _{1}$ , density  $\rho _{1}$ , shear wave velocity  $\beta _{1}=\sqrt{\mu _{1}/\rho _{1}}$ .
*   **Layer 2:** Shear modulus  $\mu _{2}$ , density  $\rho _{2}$ , shear wave velocity  $\beta _{2}=\sqrt{\mu _{2}/\rho _{2}}$ .

For Love waves to exist, the seismic energy must be trapped in the upper layer by total internal reflection. This physically requires that the shear velocity in the lower medium is faster:  $\beta _{2}>\beta _{1}$ . The Love wave phase velocity  $c$  will fall between these two bounds:  $\beta _{1}<c<\beta _{2}$ .

Since Love waves are purely horizontally polarized shear waves (SH-waves), the displacements in the  $x$  and  $z$  directions are zero. We only have displacement in the  $y$ \-direction,  $v\left(x,z,t\right)$ .

The governing wave equations for the two layers are simply the 1D scalar wave equations for shear:

$$
\nabla ^{2}v_{1}=\frac{1}{\beta _{12}}\frac{\partial ^{2}v_{1}}{\partial t^{2}}(for 0\le z\le H)
$$

$$
\nabla ^{2}v_{2}=\frac{1}{\beta _{22}}\frac{\partial ^{2}v_{2}}{\partial t^{2}}(for z>H)
$$

### 2\. The Harmonic Ansatz

We assume a harmonic wave propagating in the  $x$ \-direction with phase velocity  $c=\omega /k$ . We separate the variables to isolate the depth dependency  $V\left(z\right)$ :

$$
v_{1}\left(x,z,t\right)=V_{1}\left(z\right)e^{i\left(kx-\omega t\right)}
$$

$$
v_{2}\left(x,z,t\right)=V_{2}\left(z\right)e^{i\left(kx-\omega t\right)}
$$

Substituting these into our wave equations yields two ordinary differential equations for the depth profile. Solving them gives us the form of the wave in both layers:

**In Layer 1 (The Waveguide):** Because  $c>\beta _{1}$ , the solution is oscillatory (a standing wave in the  $z$ \-direction):

$$
V_{1}\left(z\right)=A\sin \left(zk\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)+B\cos \left(zk\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)
$$

**In Layer 2 (The Half-Space):** Because  $c<\beta _{2}$ , the wave cannot propagate downwards and must be an evanescent wave that decays exponentially with depth (where  $D=0$  to prevent infinite displacement as  $z\to \infty$ ):

$$
V_{2}\left(z\right)=Ce^{-zk\sqrt{1-\frac{c^{2}}{\beta _{22}}}}
$$

### 3\. Applying the Boundary Conditions

Now we apply the physical constraints of the architecture of the Earth:

**Condition 1: Stress-free surface at  $z=0$ ** The shear stress must vanish at the surface:  $\mu _{1}\frac{\partial V_{1}}{\partial z}=0$ . Applying this to our  $V_{1}\left(z\right)$  equation forces the sine term's coefficient to be zero ( $A=0$ ). We are left with only the cosine term:

$$
V_{1}\left(z\right)=B\cos \left(zk\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)
$$

**Condition 2: Continuity of displacement at the interface  $z=H$ ** The two layers cannot tear apart; they must move together:  $V_{1}\left(H\right)=V_{2}\left(H\right)$ .

$$
B\cos \left(Hk\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)=Ce^{-Hk\sqrt{1-\frac{c^{2}}{\beta _{22}}}}
$$

**Condition 3: Continuity of shear stress at the interface  $z=H$ ** Newton's third law requires the shear stress across the boundary to be equal:  $\mu _{1}\frac{\partial V_{1}}{\partial z}=\mu _{2}\frac{\partial V_{2}}{\partial z}$ .

$$
-\mu _{1}B\left(k\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)\sin \left(Hk\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)=-\mu _{2}C\left(k\sqrt{1-\frac{c^{2}}{\beta _{22}}}\right)e^{-Hk\sqrt{1-\frac{c^{2}}{\beta _{22}}}}
$$

### 4\. The Love Wave Dispersion Relation

To find the solution, we divide the stress continuity equation (Condition 3) by the displacement continuity equation (Condition 2). The amplitudes  $B$  and  $C$  cancel out, and the  $\sin /\cos$  becomes a tangent, leaving us with the exact analytical solution for Love waves:

$$
\tan \left(kH\sqrt{\frac{c^{2}}{\beta _{12}}-1}\right)=\frac{\mu _{2}\sqrt{1-\frac{c^{2}}{\beta _{22}}}}{\mu _{1}\sqrt{\frac{c^{2}}{\beta _{12}}-1}}
$$

### The Physical Consequence for Structures: Dispersion

Unlike the Rayleigh wave equation—which gave us a single constant velocity independent of frequency—the Love wave equation is a transcendental equation where the phase velocity  $c$  is inextricably linked to the wavenumber  $k$  (and thus the frequency  $\omega$ ).

This means **Love waves are dispersive**. Low-frequency (long wavelength) Love waves penetrate deeper into the faster mantle and travel faster. High-frequency (short wavelength) waves are trapped in the slower crust and travel slower.

In structural terms, a single seismic impulse from a fault line spreads out over time as it travels. Instead of a building experiencing a single sharp jolt, the dispersion stretches the energy into a prolonged, rolling earthquake sequence of varying frequencies. This greatly increases the probability of the ground motion sweeping through the fundamental resonant frequency of a structure.

Here is a script to visualize this dispersion. It computes the roots of the transcendental equation to plot the fundamental mode of the Love wave velocity curve.

```
#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Toggle language: 'EN' or 'ES'
LANG = 'EN'

#%%
# Define geological parameters
beta1 = 3000.0  # Crustal shear velocity (m/s)
beta2 = 4500.0  # Mantle shear velocity (m/s)
mu1 = 3.0e10    # Crustal shear modulus (Pa)
mu2 = 6.0e10    # Mantle shear modulus (Pa)

# Array of non-dimensional frequency (k*H)
kH_vals = np.linspace(0.1, 15, 200)
love_velocities = []

# The transcendental Love wave function to find roots for
def love_eq(c, kH, b1, b2, m1, m2):
    left_side = np.tan(kH * np.sqrt((c**2 / b1**2) - 1))
    right_side = (m2 * np.sqrt(1 - (c**2 / b2**2))) / (m1 * np.sqrt((c**2 / b1**2) - 1))
    return left_side - right_side

# Solve for phase velocity c at each frequency
for kH in kH_vals:
    try:
        # Search for a root between the crustal and mantle velocities
        # We add a small epsilon to avoid division by zero at the boundaries
        c_root = brentq(love_eq, beta1 + 1e-5, beta2 - 1e-5, args=(kH, beta1, beta2, mu1, mu2))
        love_velocities.append(c_root)
    except ValueError:
        love_velocities.append(np.nan)

#%%
# Bilingual dictionaries
text = {
    'EN': {
        'title': "Love Wave Dispersion Curve (Fundamental Mode)",
        'xlabel': "Non-dimensional Frequency ($kH$)",
        'ylabel': "Phase Velocity $c$ (m/s)",
        'beta1_label': "Crustal Velocity ($\\beta_1$)",
        'beta2_label': "Mantle Velocity ($\\beta_2$)",
        'curve_label': "Love Wave Velocity $c(kH)$"
    },
    'ES': {
        'title': "Curva de Dispersión de la Onda de Love (Modo Fundamental)",
        'xlabel': "Frecuencia Adimensional ($kH$)",
        'ylabel': "Velocidad de Fase $c$ (m/s)",
        'beta1_label': "Velocidad de la Corteza ($\\beta_1$)",
        'beta2_label': "Velocidad del Manto ($\\beta_2$)",
        'curve_label': "Velocidad de la Onda de Love $c(kH)$"
    }
}

#%%
# Plotting
plt.figure(figsize=(9, 6))
plt.plot(kH_vals, love_velocities, color='teal', linewidth=2.5, label=text[LANG]['curve_label'])

# Asymptotes representing the bounds
plt.axhline(y=beta2, color='darkred', linestyle='--', alpha=0.7, label=text[LANG]['beta2_label'])
plt.axhline(y=beta1, color='goldenrod', linestyle='--', alpha=0.7, label=text[LANG]['beta1_label'])

plt.title(text[LANG]['title'], fontsize=14, fontweight='bold')
plt.xlabel(text[LANG]['xlabel'], fontsize=12)
plt.ylabel(text[LANG]['ylabel'], fontsize=12)
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc='upper right')
plt.xlim(0, 15)
plt.ylim(2800, 4700)
plt.tight_layout()
plt.show()
#%%
```

Notice how at very high frequencies (large  $kH$ ), the energy is tightly confined to the upper layer, and the velocity approaches the crustal velocity  $\beta _{1}$ . At low frequencies, the wave samples deep into the mantle, and the velocity approaches  $\beta _{2}$ .

Would you like to look at the process for integrating these analytical limits—both the Rayleigh roots and the Love dispersion relations—into the boundary conditions of a numerical FDTD solver?

## Prompt:
Is there a good quality visualization of these 4 waves already made, maybe in openFoam?

## Response:
Gemini said
-----------

While OpenFOAM is an incredibly powerful open-source tool, it is fundamentally a Computational Fluid Dynamics (CFD) solver optimized for the Navier-Stokes equations using the Finite Volume Method. While it does have extensions for solid mechanics and Fluid-Structure Interaction (such as the `FOAMySees` tool, which is used to simulate post-earthquake tsunami impacts on buildings), it is rarely the tool of choice to visualize the pure elastodynamic propagation of P, S, Rayleigh, and Love waves.

### The Open-Source Standards for Elastodynamics

To get high-fidelity, 3D visualizations of the Navier-Cauchy equations, researchers rely on specialized open-source wave propagation solvers that utilize numerical methods specifically optimized for elastodynamics:

*   **SPECFEM3D:** This is arguably the gold standard for 3D seismic wave visualization. It solves the wave equations using the Spectral Element Method (SEM) on complex hexahedral meshes. It beautifully captures the exact phenomena we discussed: how volumetric P and S waves perfectly transition into dispersive Love and Rayleigh surface waves when encountering topographical variations, sedimentary basins, and structural foundations.
*   **SeisSol:** This is another state-of-the-art open-source earthquake simulator. It uses the Arbitrary High-Order Derivative Discontinuous Galerkin (ADER-DG) method on unstructured tetrahedral meshes. It excels at visualizing dynamic fault ruptures and the highly complex 3D wavefields that subsequently interact with the built environment.

### Where to Find Ready-Made Visualizations

If you want to immediately view high-quality visualizations of these four distinct waves without compiling a massive C++ or Fortran solver:

*   **IRIS (Incorporated Research Institutions for Seismology):** They host a suite of educational tools (like the _swaves_ app) and pre-rendered 3D animations that perfectly isolate the distinct particle motions (longitudinal, transverse, rolling, and shearing) of all four waves interacting with a free surface.
*   **Scientific Repositories:** The documentation pages and publication galleries for both SPECFEM3D and SeisSol feature extensive video libraries of large-scale simulations (like the 2019 Ridgecrest earthquake) where the distinct separation of body and surface wavefronts is visually stunning.

Would you like me to write a custom Python script utilizing a 3D visualization library like PyVista to generate an interactive, bilingual (EN/ES) 3D animation of these surface wave particle motions for you to run locally?

* * *

[High Performance Seismic Wave Propagation with SPECFEM3D](https://www.youtube.com/watch?v=pfT2xOLz_Vc)
 This presentation showcases how the SPECFEM3D package utilizes the spectral element method to generate high-fidelity seismic wave visualizations on large arbitrary domains.

> [High Performance Seismic Wave Propagation with SPECFEM3D - YouTube - insideHPC Report · 3.2K views](https://www.youtube.com/watch?v=pfT2xOLz_Vc)



---
Powered by [Gemini Exporter](https://www.ai-chat-exporter.com)