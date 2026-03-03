#= 
-------------------------------------------------------------------------
Example 3: Solid Sphere in Own Gravity (Landau & Lifshitz, §7 Prob 3)
-------------------------------------------------------------------------
=#

"""
    gravitating_sphere(r, params)

Calculates the radial displacement and central pressure of a self-gravitating sphere.

# Arguments
- `r`: Radial coordinate [m]
- `params`: NamedTuple with:
    - `R`:   Sphere radius [m]
    - `rho`: Density [kg/m³]
    - `g`:   Surface gravity [m/s²]
    - `E`:   Young's modulus [Pa]
    - `nu`:  Poisson's ratio [-] (sigma in text)

# Returns
- `u_r`: Radial displacement at r [m]
- `p_center`: Pressure at the center of the sphere [Pa]
"""
function gravitating_sphere(r::Float64, params)
    R, ρ, g = params.R, params.rho, params.g
    E, σ    = params.E, params.nu

    # Pre-calculate material factor
    # Factor = (g rho R (1-2σ)(1+σ)) / (10 E (1-σ))
    mat_factor = (g * ρ * R * (1 - 2*σ) * (1 + σ)) / (10 * E * (1 - σ))

    # Geometric term: (3-σ)/(1+σ) - (r/R)^2
    geom_term = ((3 - σ) / (1 + σ)) - (r / R)^2

    # Displacement u_r
    u_r = -mat_factor * r * geom_term

    # Pressure at center (compressive stress)
    # p = (3-σ)/(10(1-σ)) * g rho R
    p_center = ((3 - σ) / (10 * (1 - σ))) * g * ρ * R

    return u_r, p_center
end

# --- Test Case: Earth-like Planetoid (Iron) ---
planet_params = (
    R   = 1000e3,   # [m] (1000 km radius)
    rho = 7800.0,   # [kg/m³]
    g   = 1.0,      # [m/s²] (Assumed surface gravity)
    E   = 2.1e11,   # [Pa]
    nu  = 0.3       # [-]
)

u_surf, p_core = gravitating_sphere(planet_params.R, planet_params)

println("--- Example 3: Self-Gravitating Sphere ---")
println("Surface Displacement: $(round(u_surf, digits=4)) [m] (Negative means shrinkage)")
println("Core Pressure:        $(round(p_core / 1e9, digits=2)) [GPa]")