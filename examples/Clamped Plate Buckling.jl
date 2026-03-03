# import Pkg; Pkg.add("Arpack"); Pkg.add("SparseArrays")
using LinearAlgebra, SparseArrays, Arpack, Plots

#= 
-------------------------------------------------------------------------
 Numerical Plate Buckling (Finite Difference Method)
 Case: Clamped Edges (Fixed Displacement, Fixed Rotation)
 Solving: D * grad^4(w) = -N_x * d^2w/dx^2
-------------------------------------------------------------------------
=#

function solve_plate_buckling(Lx, Ly, Nx, Ny)
    # Grid steps
    hx = Lx / (Nx + 1)
    hy = Ly / (Ny + 1)

    # --- 1. Construct 1D Operators ---
    # We solve for INTERNAL nodes only.
    # For Clamped BCs (u=0, u'=0 at boundary), the stencil near the edge changes.
    # Standard Centered Difference for 2nd deriv: [1, -2, 1]
    # Standard Centered Difference for 4th deriv: [1, -4, 6, -4, 1]
    
    # 2nd Derivative Matrix (1D)
    D2x = spdiagm(-1 => ones(Nx-1), 0 => -2*ones(Nx), 1 => ones(Nx-1)) / hx^2
    D2y = spdiagm(-1 => ones(Ny-1), 0 => -2*ones(Ny), 1 => ones(Ny-1)) / hy^2

    # 4th Derivative Matrix (1D) for CLAMPED ends
    # The boundary condition u=0 is handled by excluding boundary nodes.
    # The boundary condition u'=0 modifies the near-boundary stencil.
    # (Ghost point u_-1 = u_1) -> Stencil at i=1 becomes [7, -4, 1] instead of [1, -4, 6...]
    
    # Helper to build 4th deriv matrix
    function build_D4_clamped(N, h)
        main = 6 * ones(N)
        off1 = -4 * ones(N-1)
        off2 = 1 * ones(N-2)
        
        # Modify for Clamped BCs (u=0, u'=0)
        # Using ghost point method: u_0 = 0, u_-1 = u_1
        # Row 1: u_-1 - 4u_0 + 6u_1 - 4u_2 + u_3 -> (u_1) - 0 + 6u_1 - 4u_2 + u_3 = 7u_1 - 4u_2 + u_3
        main[1] += 1
        main[end] += 1
        
        return spdiagm(-2 => off2, -1 => off1, 0 => main, 1 => off1, 2 => off2) / h^4
    end

    D4x = build_D4_clamped(Nx, hx)
    D4y = build_D4_clamped(Ny, hy)
    
    # Identity Matrices
    Ix = sparse(LinearAlgebra.I, Nx, Nx)
    Iy = sparse(LinearAlgebra.I, Ny, Ny)

    # --- 2. Construct 2D Operators (Kronecker Products) ---
    # Biharmonic Operator: D^4 = D4x⊗Iy + 2*(D2x⊗D2y) + Ix⊗D4y
    # Note: The cross term 2*d4/dx2dy2 corresponds to D2x ⊗ D2y
    
    L_biharm = kron(D4x, Iy) + 2 * kron(D2x, D2y) + kron(Ix, D4y)
    
    # Load Operator (Compression in X): -d^2/dx^2
    # Operator is applied to w. N_x * d^2w/dx^2 -> Matrix is D2x ⊗ Iy
    L_load = -kron(D2x, Iy)

    # --- 3. Solve Eigenvalue Problem ---
    # A * u = lambda * B * u
    # We look for the SMALLEST POSITIVE eigenvalue (Lowest Critical Load)
    
    # eigs solves A*x = lambda*B*x. We want smallest magnitude.
    # sigma=0.0 finds eigenvalues near 0.
    vals, vecs = eigs(L_biharm, L_load; nev=1, which=:SM, sigma=0.0)
    
    critical_load_factor = vals[1] # This is N_x / D
    mode_shape = reshape(real(vecs[:,1]), Nx, Ny)
    
    return critical_load_factor, mode_shape
end

# Parameters
L_plate = 2.0  # Length (a)
W_plate = 1.0  # Width (b) - Aspect Ratio = 2.0
Grid_N = 40    # Resolution

# Theoretical comparison for Clamped Plate (approximate)
# For Aspect Ratio 2.0, k is approx 7.0 (higher than SS case of 4.0)
N_cr_normalized, shape = solve_plate_buckling(L_plate, W_plate, Grid_N*2, Grid_N)

println("--- Numerical Solution (Clamped Plate) ---")
println("Computed Eigenvalue (N_cr / D): $(round(N_cr_normalized, digits=4))")
# Convert back to Buckling Coefficient k for comparison: 
# sigma = k * pi^2 * D / (b^2 * h) -> N_cr = k * pi^2 * D / b^2
# k = N_cr_normalized * b^2 / pi^2
k_computed = N_cr_normalized * W_plate^2 / pi^2
println("Computed Buckling Coefficient k: $(round(k_computed, digits=2))")
println("(Note: For Clamped-Clamped a/b=2, k should be approx 6.9 - 7.2)")

# Plotting
heatmap(shape', title="Buckling Mode Shape (Clamped)", 
        aspect_ratio=:equal, c=:viridis, 
        xlabel="Length (x)", ylabel="Width (y)")