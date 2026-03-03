using Plots

#= 
-------------------------------------------------------------------------
 Plate Buckling Analysis (Box Girder Bridge)
 Calculating the Buckling Coefficient k vs Aspect Ratio (a/b)
-------------------------------------------------------------------------
=#

function buckling_coefficient(aspect_ratio, m)
    # k = (m * (b/a) + (a/b) / m)^2
    # Let x = a/b (Aspect Ratio)
    # k = (m/x + x/m)^2
    return (m ./ aspect_ratio .+ aspect_ratio ./ m).^2
end

# 1. Generate the Garland Curves
aspect_ratios = range(0.1, 5.0, length=200)
k_m1 = buckling_coefficient.(aspect_ratios, 1)
k_m2 = buckling_coefficient.(aspect_ratios, 2)
k_m3 = buckling_coefficient.(aspect_ratios, 3)
k_m4 = buckling_coefficient.(aspect_ratios, 4)

# 2. Determine the critical k (minimum of all modes)
k_crit = min.(k_m1, k_m2, k_m3, k_m4)

# 3. Calculate Critical Stress for a specific Steel Plate
# Parameters: Steel Bridge Deck
E_steel = 210e9  # Pa
nu_steel = 0.3
b_width = 1.0    # m (spacing between webs or stiffeners)
h_thick = 0.015  # m (15 mm thickness)

# Formula: sigma_cr = (k * pi^2 * E) / (12 * (1-nu^2)) * (h/b)^2
# Let's verify for a long plate (a >> b), where k_min = 4.0
k_min_long = 4.0
sigma_cr_val = (k_min_long * pi^2 * E_steel) / (12 * (1 - nu_steel^2)) * (h_thick / b_width)^2

println("--- Plate Buckling Results ---")
println("Critical Stress for 15mm plate, 1m wide: $(round(sigma_cr_val/1e6, digits=2)) MPa")
if sigma_cr_val < 350e6 # Typical Yield Strength of structural steel
    println("Status: Elastic Buckling governs (Fails before Yielding)")
else
    println("Status: Yielding governs (Plate is stocky enough)")
end

# 4. Plot
plot(aspect_ratios, [k_m1, k_m2, k_m3, k_m4], 
     label=["m=1" "m=2" "m=3" "m=4"], 
     ylim=(0, 10), lw=2, linestyle=:dash)

plot!(aspect_ratios, k_crit, 
      label="Critical k", color=:black, lw=3,
      title="Plate Buckling Coefficient k (Garland Curves)",
      xlabel="Aspect Ratio (a/b)", ylabel="Buckling Coefficient k")