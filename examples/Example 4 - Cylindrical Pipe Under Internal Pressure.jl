#= 
-------------------------------------------------------------------------
Example 4: Cylindrical Pipe (Landau & Lifshitz, §7 Prob 4)
-------------------------------------------------------------------------
=#

"""
    pipe_deformation(r, params)

Calculates displacement and hoop stress for a pipe under internal pressure.

# Arguments
- `r`: Radial coordinate [m]
- `params`: NamedTuple with:
    - `R1`, `R2`: Inner/Outer radii [m]
    - `p`:        Internal pressure [Pa]
    - `E`:        Young's modulus [Pa]
    - `nu`:       Poisson's ratio [-]

# Returns
- `u`:        Radial displacement [m]
- `sigma_phi`: Hoop stress [Pa]
"""
function pipe_deformation(r::Float64, params)
    R1, R2 = params.R1, params.R2
    p, E, σ = params.p, params.E, params.nu

    if r < R1 || r > R2
        error("r out of bounds")
    end

    # Common factor: P * R1^2 / (R2^2 - R1^2)
    geom_factor = (p * R1^2) / (R2^2 - R1^2)

    # Displacement u
    # u = ((1+σ)/E) * geom_factor * [ (1-2σ)r + R2^2/r ]
    term_bracket = (1 - 2*σ)*r + (R2^2 / r)
    u = ((1 + σ) / E) * geom_factor * term_bracket

    # Stress sigma_phi (Hoop stress)
    # Formula derived from text: sigma_phi = geom_factor * (1 + R2^2/r^2)
    sigma_phi = geom_factor * (1 + (R2^2 / r^2))

    return u, sigma_phi
end