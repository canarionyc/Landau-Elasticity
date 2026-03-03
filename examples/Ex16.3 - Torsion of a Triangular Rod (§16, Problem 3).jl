#= 
-------------------------------------------------------------------------
 §16 Prob 3: Torsion of Equilateral Triangle
-------------------------------------------------------------------------
=#

function torsion_triangle(a, mu, grid_size=100)
    # Torsional Rigidity
    C = (sqrt(3) * mu * a^4) / 80
    
    # Grid setup (Triangle height h = a*sqrt(3)/2)
    h = a * sqrt(3) / 2
    # Centroid is at h/3 from base. Let's center coordinates there.
    # Vertices: (-a/2, -h/3), (a/2, -h/3), (0, 2h/3)
    
    xs = range(-a/2, a/2, length=grid_size)
    ys = range(-h/3, 2*h/3, length=grid_size)
    
    # Stress vector magnitude |tau|
    stress_mag = zeros(grid_size, grid_size)
    
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys)
        # Check if inside triangle (barycentric or plane equations)
        # Using simple line equations for boundaries
        in_tri = (y >= -h/3) && (y <= -sqrt(3)*x + 2*h/3) && (y <= sqrt(3)*x + 2*h/3)
        
        if in_tri
            # Derivatives of psi = y(3x^2 - y^2) / 6a
            dpsi_dx = (y * 6 * x) / (6 * a)
            dpsi_dy = (3 * x^2 - 3 * y^2) / (6 * a)
            
            sx = mu * (dpsi_dx - y)
            sy = mu * (dpsi_dy + x)
            stress_mag[j, i] = sqrt(sx^2 + sy^2)
        else
            stress_mag[j, i] = NaN
        end
    end
    
    return C, xs, ys, stress_mag
end

C_tri, x_tri, y_tri, stress_field = torsion_triangle(0.1, 80e9)

heatmap(x_tri, y_tri, stress_field, title="Shear Stress Magnitude (Triangular Rod)",
        aspect_ratio=:equal, c=:inferno)