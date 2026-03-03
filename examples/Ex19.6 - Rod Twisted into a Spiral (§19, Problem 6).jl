using LinearAlgebra

#= 
-------------------------------------------------------------------------
 §19 Prob 6: Rod Twisted into a Spiral
-------------------------------------------------------------------------
=#

function spiral_rod_forces(R, alpha_deg, E, mu, r_rod)
    alpha = deg2rad(alpha_deg)
    
    # Section properties (Circular rod radius r_rod)
    I = (pi * r_rod^4) / 4
    J = (pi * r_rod^4) / 2 # Polar moment
    C = mu * J             # Torsional rigidity
    EI = E * I             # Flexural rigidity
    
    # Torsion of the spiral geometry
    # The twist 'tau' here refers to the geometric torsion of the curve
    tau_geom = (sin(alpha) * cos(alpha)) / R
    
    # 1. Calculate Forces and Moments required to hold this shape
    # Force F (along Z axis)
    term1 = (C * tau_geom / R) * sin(alpha)
    term2 = (EI / R^2) * (cos(alpha)^2) * sin(alpha)
    F_z = term1 - term2
    
    # Moment M_z
    M_z = C * tau_geom * sin(alpha) + (EI / R) * cos(alpha)^3
    
    # 2. Generate 3D Geometry for plotting
    phi = range(0, 4*pi, length=200) # 2 turns
    x_vals = R .* cos.(phi)
    y_vals = R .* sin.(phi)
    z_vals = (R * tan(alpha)) .* phi
    
    return F_z, M_z, x_vals, y_vals, z_vals
end

# Parameters
R_cyl = 1.0      # Helix radius
alpha = 30.0     # Pitch angle
radius_rod = 0.05
E_steel = 210e9
G_steel = 80e9

F_req, M_req, xs, ys, zs = spiral_rod_forces(R_cyl, alpha, E_steel, G_steel, radius_rod)

println("Force Fz required: $(round(F_req, digits=2)) N")
println("Moment Mz required: $(round(M_req, digits=2)) N*m")

# 3D Plot
plot(xs, ys, zs, label="Deformed Rod", lw=3, title="Rod Twisted into Spiral",
     xlabel="x", ylabel="y", zlabel="z", aspect_ratio=:equal)

# Visualizing the "Spring Constant" from the spiral rod equations
alphas = range(5, 60, length=50) # Varying the pitch
forces = [spiral_rod_forces(1.0, a, 210e9, 80e9, 0.01)[1] for a in alphas]

plot(alphas, forces, title="Force vs. Pitch Angle (Spring Effect)",
     xlabel="Pitch Angle alpha (degrees)", ylabel="Required Axial Force (N)",
     lw=2, legend=false)