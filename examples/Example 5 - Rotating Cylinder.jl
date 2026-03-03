#= 
-------------------------------------------------------------------------
Example 5: Rotating Cylinder (Landau & Lifshitz, §7 Prob 5)
-------------------------------------------------------------------------
=#

function rotating_cylinder(r::Float64, params)
    R, ρ, Ω = params.R, params.rho, params.Omega
    E, σ    = params.E, params.nu

    # Pre-factor
    # coeff = (rho Omega^2 (1+σ)(1-2σ)) / (8 E (1-σ))
    coeff = (ρ * Ω^2 * (1 + σ) * (1 - 2*σ)) / (8 * E * (1 - σ))

    # Displacement
    u = coeff * r * ((3 - 2*σ)*R^2 - r^2)

    return u
end