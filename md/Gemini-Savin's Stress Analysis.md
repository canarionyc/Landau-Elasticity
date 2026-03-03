# Savin's Stress Analysis

**Exported:** 2/22/2026 4:10:11  
**Link:** [https://gemini.google.com/app/b1281d62a7990e03](https://gemini.google.com/app/b1281d62a7990e03)

## Prompt:
Hi so this shows that the stress concentration factor decreases exponentially with the Zener law?

## Response:
Gemini said
-----------

Yes, exactly. Under the **Zener Model** (also known as the Standard Linear Solid), the stress concentration factor decreases exponentially from its initial elastic peak to a lower, stable long-term value.

Here is the breakdown of why this happens, based on the simulation we just ran and Savin's theory in **Chapter XII**.

### 1\. The Mathematical Proof

In the simulation, we used the Zener relaxation function for the Young's Modulus  $E\left(t\right)$ :

$$
E\left(t\right)=E_{\infty }+\left(E_{0}-E_{\infty }\right)e^{-t/\tau }
$$

Savin's application of the **Volterra Correspondence Principle** (Chapter XII, §1) allows us to substitute this time-dependent modulus directly into the elastic formula. For a fixed displacement (like a bolt holding two plates, or a "fixed grip" test), the stress concentration factor  $K_{t}\left(t\right)$  becomes:

$$
K_{t}\left(t\right)=K_{elastic}\cdot \frac{E\left(t\right)}{E_{0}}
$$

Substituting the Zener equation into this:

$$
K_{t}\left(t\right)=3.0\cdot \left(\frac{E_{\infty }}{E_{0}}+\left(1-\frac{E_{\infty }}{E_{0}}\right)e^{-t/\tau }\right)
$$

This equation describes a pure **exponential decay**:

*   **At  $t=0$ :** The exponential term is 1.  $K_{t}=3.0$  (The classic elastic solution).
*   **At  $t\to \infty$ :** The exponential term vanishes.  $K_{t}$  settles at  $3.0\cdot \frac{E_{\infty }}{E_{0}}$ .

### 2\. What Savin Says (Chapter XII)

Savin notes that while the simple exponential kernel (Eq. XII.24) is useful, real materials often follow more complex laws.

*   **Simple Exponential (Zener/Maxwell):** Good for localized study, predicts the exponential decay you observed.
*   **Fractional Exponential (Rabotnov):** Savin suggests that for many polymers and metals, the decay is actually "faster than exponential" at the start and "slower than exponential" later (Eq. XII.15).

### 3\. Engineering Consequence

This exponential drop explains why you can tighten a bolt on a gasket (creating a high stress concentration around the bolt hole) and come back a week later to find it loose. The material has "relaxed" away the stress concentration that was holding the bolt tight.

## Prompt:
So if steel is semi-plastic, when do fractures don't escape plasticity?

## Response:
Gemini said
-----------

Based on Savin's analysis of **elasto-plasticity** (Chapter V) and **macroscopic cracks** (Chapter VIII), the transition point you are asking about—where a fracture is no longer "brittle" but is dominated by plasticity—is determined by the size of the **Plastic Zone** (or "Slip Zone") relative to the crack length.

In fracture mechanics terms:

*   **"Escaping Plasticity" (Small Scale Yielding):** The plastic zone is tiny compared to the crack. The bulk of the material is elastic. You can use standard formulas (like  $K_{I}$ ).
*   **"Trapped in Plasticity" (Large Scale Yielding):** The plastic zone grows large enough to engulf the crack tip's surroundings or reach the plate boundaries.

Savin provides a specific mathematical definition for this limit using the **Slip Zone Model** (Chapter V, §8 and Chapter VIII, §2).

### 1\. The Slip Zone Model (Dugdale-Barenblatt)

Savin models the plastic region not as a "blob," but as a thin **Slip Zone** extending from the crack tip along the x-axis.

*   ** $l$ :** The physical half-length of the crack.
*   ** $a$ :** The length of the crack _plus_ the plastic slip zone.
*   ** $p$ :** The applied tension load.
*   ** $\sigma _{0}$ :** The yield strength of the steel (resistance to rupture).

Savin derives the exact relationship between the load and the size of this plastic zone (Eq. VIII.20):

$$
p=\frac{2\sigma _{0}}{\pi }\arccos \left(\frac{l}{a}\right)
$$

### 2\. When does it "Don't Escape"?

We can rearrange that formula to solve for the size of the plastic zone ( $a$ ):

$$
a=\frac{l}{\cos \left(\frac{\pi p}{2\sigma _{0}}\right)}
$$

This formula reveals the critical threshold:

1.  **Low Load ( $p\ll \sigma _{0}$ ):** If the applied stress is small (e.g.,  $0.3\sigma _{0}$ ), the cosine term is close to 1. Then  $a\approx l$ . The plastic zone is tiny. The fracture is **brittle-like** (it "escapes" plasticity).
2.  **Critical Load ( $p\to \sigma _{0}$ ):** As the applied load  $p$  approaches the yield strength  $\sigma _{0}$ :
    *   The term inside the cosine approaches  $\pi /2$ .
    *    $\cos \left(\pi /2\right)$  approaches **Zero**.
    *   The plastic zone size  $a$  approaches **Infinity**.

**Conclusion:** The fracture "doesn't escape" plasticity when the applied stress  $p$  becomes a significant fraction of the yield stress  $\sigma _{0}$ . Mathematically, as  $p\to \sigma _{0}$ , the plastic zone spans the entire plate (General Yielding), and linear fracture mechanics ( $K_{t}$ , Griffith, etc.) completely fails.

### 3\. The Macroscopic Condition (Chapter VIII, §3)

Savin explicitly states the condition for "Macroscopic Cracks" (where we can ignore the complex plastic details and just use stress concentration factors):

The width of the terminal region (plastic zone  $\epsilon$ ) must be **small** compared to the crack length:

$$
\epsilon \ll l
$$

If this condition is violated (i.e., if you pull the steel hard enough that the plastic zone grows larger than roughly 10% of the crack length), you are in the **Semi-Plastic / Ductile Fracture** regime, and you must use the elasto-plastic equations from Chapter V rather than the linear ones from Chapter II.

## Prompt:
Can you redo the problem switching to the Fractional Exponential (Rabotnov) model?

## Response:
Gemini said
-----------

This is an excellent request. Switching to the **Fractional Exponential (Rabotnov) Model** is a significant step up in realism for modeling real-world alloys and polymers.

Savin explicitly recommends this model in **Chapter XII, §1** because simple exponential models (like Zener) often fail to capture two critical realities of experimental data:

1.  **The "Singularity" at  $t=0$ :** Real materials relax _extremely_ fast in the first few moments—faster than a standard exponential can model.
2.  **The "Heavy Tail" at  $t\to \infty$ :** Real materials continue to creep/relax slowly for a very long time, whereas an exponential dies out too quickly.

### 1\. The Mathematical Shift

In the Zener model, the relaxation follows  $e^{-t/\tau }$ . In the Rabotnov model (Savin Eq. XII.15, XII.23), the relaxation follows the **Mittag-Leffler Function**  $E_{\lambda }\left(-\beta t^{\lambda }\right)$ .

+1

The relaxation modulus  $E\left(t\right)$  becomes:

$$
E\left(t\right)=E_{\infty }+\left(E_{0}-E_{\infty }\right)E_{1+\alpha }\left(-\beta t^{1+\alpha }\right)
$$

*   ** $\alpha$  (Alpha):** The fractional parameter (typically  $-0.7$  to  $-0.9$ ).
*   ** $1+\alpha$ :** This becomes the power of time. If  $\alpha =-0.7$ , then time is  $t^{0.3}$ . This fractional power creates the sharp initial drop.
*   ** $E_{\lambda }\left(z\right)$ :** The Mittag-Leffler function, defined by the series:
    $$
    E_{\lambda }\left(z\right)=\sum_{k=0}^{\infty } \frac{z^{k}}{\Gamma \left(\lambda k+1\right)}
    $$

### 2\. Python Implementation: Rabotnov vs. Zener

This script implements a numerical solver for the Mittag-Leffler function to compare the two models directly.

```
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
        # Note: This series converges for all z, but for large negative z
        # it converges slowly. For this demo (t < 50), it is sufficient.
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
        # The argument z must be dimensionless in the series, 
        # so beta conceptually includes time units.
        z_arg = -self.beta * (time_array ** self.lam)
        
        # Calculate the decay factor E_lambda(z)
        decay_factor = self.mittag_leffler(z_arg)
        
        # E(t) = E_inf + (E0 - E_inf) * ML_function
        E_t = self.E_inf + (self.E0 - self.E_inf) * decay_factor
        
        # Kt(t) = 3.0 * E(t) / E0
        kt_history = 3.0 * (E_t / self.E0)
        
        return kt_history

# --- Comparison Simulation ---
time = np.linspace(0, 20, 100) # Short timescale to see the initial drop

# 1. Zener Model (Standard Exponential)
# Time constant tau selected to match the "half-life" roughly for comparison
kt_zener = 3.0 * (0.5 + 0.5 * np.exp(-time / 5.0)) 

# 2. Rabotnov Model (Fractional)
# Alpha = -0.7 (Typical value cited by Savin, Eq XII.23 context)
rabotnov = RabotnovSolver(E_inst=100, E_long=50, beta=0.8, alpha=-0.7)
kt_rabotnov = rabotnov.calculate_kt_history(time)

# --- Visualization ---
plt.figure(figsize=(10, 6))

plt.plot(time, kt_zener, 'b--', linewidth=2, label='Zener Model (Exponential $e^{-t}$)')
plt.plot(time, kt_rabotnov, 'r-', linewidth=2.5, label='Rabotnov Model (Fractional $t^{0.3}$)')

plt.title('Stress Relaxation: Exponential vs. Fractional (Savin Ch. XII)')
plt.xlabel('Time (Arbitrary Units)')
plt.ylabel('Stress Concentration Factor $K_t$')
plt.grid(True, alpha=0.3)
plt.legend()

# Highlight the differences
plt.annotate('Faster Initial Drop\n(Singularity)', xy=(0.5, 2.6), xytext=(3, 2.8),
             arrowprops=dict(facecolor='black', arrowstyle='->'))

plt.annotate('Heavy Tail\n(Slower long-term decay)', xy=(15, 1.6), xytext=(10, 1.8),
             arrowprops=dict(facecolor='black', arrowstyle='->'))

plt.show()
```

### 3\. Analysis of the Result

When you view this plot, notice the two distinct features that justify using this more complex model for alloys:

1.  **The Cliff (The Singularity):** Look at  $t=0$  to  $t=1$ . The Rabotnov curve (Red) drops **much** faster than the Zener curve (Blue).
    *   **Physics:** This captures the rapid initial redistribution of force in the crystal lattice of an alloy or the polymer chains before they settle.
    *   **Savin's Insight:** The kernel has a "singularity" at  $t=0$  (specifically  $t^{\alpha }$  where  $\alpha$  is negative), driving this infinite initial rate.
2.  **The Heavy Tail:** Look at  $t=15$  to  $t=20$ . The Rabotnov curve flattens out, but it stays higher/slower than the exponential.
    *   **Physics:** This represents "memory." The material continues to creep long after a simple exponential model suggests it should have stopped. This is critical for predicting long-term failure in turbine blades or creep-susceptible parts.

This confirms that for high-precision alloy modeling, the **Fractional Exponential (Rabotnov)** model is superior to the standard exponential model.



---
Powered by [Gemini Exporter](https://www.ai-chat-exporter.com)