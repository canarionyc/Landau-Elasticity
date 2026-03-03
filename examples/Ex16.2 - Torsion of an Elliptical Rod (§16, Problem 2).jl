using Plots

#= 
-------------------------------------------------------------------------
 §16 Prob 2: Torsion of Elliptical Rod
-------------------------------------------------------------------------
=#

function torsion_ellipse(a, b, mu, grid_size=100)
    # 1. Calculate Torsional Rigidity C analytically
    C = pi * mu * (a^3 * b^3) / (a^2 + b^2)
    
    # 2. Define the Stress Function Chi field
    xs = range(-a, a, length=grid_size)
    ys = range(-b, b, length=grid_size)
    chi = zeros(grid_size, grid_size)
    
    coeff = (a^2 * b^2) / (2 * (a^2 + b^2))
    
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
        val = 1 - (x^2/a^2) - (y^2/b^2)
        if val >= 0
            chi[j, i] = coeff * val
        else
            chi[j, i] = NaN # Outside the rod
        end
    end
    
    return C, xs, ys, chi
end

# Parameters
a, b = 2.0, 1.0 # Semi-axes
mu_val = 80e9   # Shear Modulus (Pa)

C_val, x_grid, y_grid, chi_field = torsion_ellipse(a, b, mu_val)

println("Torsional Rigidity C: $(C_val) N*m^2")

# Plotting the Stress Function (Contour Plot)
contour(x_grid, y_grid, chi_field, fill=true, title="Stress Function Chi (Elliptical Rod)", 
        xlabel="x", ylabel="y", aspect_ratio=:equal, c=:viridis)