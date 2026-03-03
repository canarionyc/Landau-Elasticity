#= 
-------------------------------------------------------------------------
Example 12: Stress Concentration (Spherical Cavity)
-------------------------------------------------------------------------
=#

"""
    cavity_stress_concentration(sigma_0, nu)

Calculates the maximum stress at the equator of a spherical cavity 
in an infinite medium under uniaxial tension sigma_0.

# Arguments
- `sigma_0`: Applied tension at infinity [Pa]
- `nu`:      Poisson's ratio [-]

# Returns
- `sigma_max`: Maximum stress at the cavity surface [Pa]
- `k`:         Stress concentration factor (sigma_max / sigma_0)
"""
function cavity_stress_concentration(sigma_0::Float64, nu::Float64)
    # Formula: k = (27 - 15σ) / (2 * (7 - 5σ))
    numerator   = 27.0 - 15.0 * nu
    denominator = 2.0 * (7.0 - 5.0 * nu)
    
    k = numerator / denominator
    sigma_max = k * sigma_0
    
    return sigma_max, k
end

# --- Test Case ---
# Steel (nu = 0.3)
s0 = 100.0e6 # 100 MPa
s_max, k_factor = cavity_stress_concentration(s0, 0.3)

println("\n--- Example 12: Spherical Cavity Stress Concentration ---")
println("Applied Stress: $(s0/1e6) MPa")
println("Max Stress:     $(round(s_max/1e6, digits=2)) MPa")
println("Concentration Factor k: $(round(k_factor, digits=3))")