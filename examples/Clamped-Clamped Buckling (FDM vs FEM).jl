using LinearAlgebra, SparseArrays, Plots, Arpack

# --- PHYSICAL PARAMETERS ---
L = 2.0          # Length
EI = 1000.0      # Flexural Rigidity
Target_P = 4 * pi^2 * EI / L^2  # Analytical Solution (approx 9869.6)

# =========================================================
# METHOD 1: FINITE DIFFERENCE (FDM)
# =========================================================
function solve_fdm(N)
    h = L / (N + 1)
    
    # 4th Derivative Stencil: [1, -4, 6, -4, 1]
    # We construct the matrix for internal nodes only
    # BCs: u=0 at ends implies boundary nodes are 0.
    # Clamped (u'=0) implies ghost points u_{-1} = u_{1}.
    
    # Construct diagonals for the Biharmonic Operator (d^4/dx^4)
    main = 6 * ones(N)
    main[1] += 1; main[end] += 1  # Clamped correction
    off1 = -4 * ones(N-1)
    off2 = 1 * ones(N-2)
    A = (spdiagm(-2=>off2, -1=>off1, 0=>main, 1=>off1, 2=>off2))

    # Construct diagonals for the Harmonic Operator (d^2/dx^2) for Load
    main_B = -2 * ones(N)
    off1_B = 1 * ones(N-1)
    B = (spdiagm(-1=>off1_B, 0=>main_B, 1=>off1_B))
    
    # Solve: (EI/h^4) * A * w = - (P/h^2) * B * w
    # Generalized Eigenvalue: A*w = lambda * (-B) * w
    # where lambda = (P * h^2) / EI
    
    vals, vecs = eigs(A, -B; nev=1, which=:SM) 
    P_crit = vals[1] * EI / h^2
    return P_crit
end

# =========================================================
# METHOD 2: FINITE ELEMENT METHOD (FEM)
# =========================================================
function solve_fem(n_elements)
    # Length of one element
    le = L / n_elements
    
    # Degrees of Freedom (DoF)
    # Each node has 2 DoFs: Vertical Displacement (w) and Rotation (theta)
    # Total DoFs = 2 * (n_elements + 1)
    n_dof = 2 * (n_elements + 1)
    
    # Global Matrices
    K = spzeros(n_dof, n_dof) # Stiffness
    Kg = spzeros(n_dof, n_dof) # Geometric Stiffness (Load)
    
    # --- Element Matrices (Euler-Bernoulli Beam) ---
    # These are standard derivation results for cubic shape functions
    
    # Stiffness Element (Bending Energy)
    # Factor: EI / le^3
    k_local = [ 12   6*le    -12   6*le;
                6*le 4*le^2  -6*le 2*le^2;
               -12  -6*le     12  -6*le;
                6*le 2*le^2  -6*le 4*le^2 ] * (EI / le^3)
                
    # Geometric Element (Work done by P)
    # Factor: 1 / (30 * le)
    kg_local = [ 36     3*le    -36     3*le;
                 3*le   4*le^2  -3*le  -le^2;
                -36    -3*le     36    -3*le;
                 3*le  -le^2    -3*le   4*le^2 ] * (1.0 / (30 * le))

    # --- Assembly ---
    for i in 1:n_elements
        # Indices for the 4 DoFs of this element in the global matrix
        # Node i and Node i+1
        # Node i DoFs: 2i-1, 2i
        idx = [2*i-1, 2*i, 2*i+1, 2*i+2]
        
        K[idx, idx] += k_local
        Kg[idx, idx] += kg_local
    end
    
    # --- Apply Boundary Conditions ---
    # Clamped at both ends: w=0, theta=0 at Node 1 and Node end
    # We constrain DoFs: 1, 2 (Left) and n_dof-1, n_dof (Right)
    fixed_dofs = [1, 2, n_dof-1, n_dof]
    free_dofs = setdiff(1:n_dof, fixed_dofs)
    
    K_reduced = K[free_dofs, free_dofs]
    Kg_reduced = Kg[free_dofs, free_dofs]
    
    # Solve K * u = P * Kg * u
    vals, vecs = eigs(K_reduced, Kg_reduced; nev=1, which=:SM)
    return vals[1]
end

# --- RUN COMPARISON ---
println("Analytical Target: $Target_P")

# Run with coarse grid (low number of points) to see which converges faster
N_coarse = 6

P_fdm = solve_fdm(N_coarse)
P_fem = solve_fem(N_coarse)

println("\n--- Comparison (N = $N_coarse elements/points) ---")
println("FDM Result: $P_fdm \t Error: $(round(100*abs(P_fdm-Target_P)/Target_P, digits=2))%")
println("FEM Result: $P_fem \t Error: $(round(100*abs(P_fem-Target_P)/Target_P, digits=2))%")

println("\nFEM is typically vastly more accurate for the same number of points in structural problems.")

# Analytical Target: 9869.604401089358


#= --- Comparison (N = 6 elements/points) ---
FDM Result: 9224.499854461037    Error: 6.54%
FEM Result: 9885.211838001385    Error: 0.16%

FEM is typically vastly more accurate for the same number of points in structural problems.

Analytical Target: 9869.604401089358

--- Comparison (N = 6 elements/points) ---
FDM Result: 9224.49985446103     Error: 6.54%
FEM Result: 9885.211838001389    Error: 0.16%

FEM is typically vastly more accurate for the same number of points in structural problems.

Analytical Target: 9869.604401089358

--- Comparison (N = 6 elements/points) ---
FDM Result: 9224.499854461043    Error: 6.54%
FEM Result: 9885.211838001398    Error: 0.16%

FEM is typically vastly more accurate for the same number of points in structural problems.
 =#