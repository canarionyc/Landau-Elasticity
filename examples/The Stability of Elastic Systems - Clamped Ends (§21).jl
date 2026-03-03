using DifferentialEquations, Plots

#= 
-------------------------------------------------------------------------
 §21: Buckling with Clamped Ends
 Boundary Conditions: X(0)=0, X'(0)=0, X(L)=0, X'(L)=0
-------------------------------------------------------------------------
=#

function buckling_clamped(L, T_guess, EI)
    # ODE: X'''' + k^2 * X'' = 0
    # Let y[1]=X, y[2]=X', y[3]=X'', y[4]=X'''
    k2 = T_guess / EI
    
    function buck_ode!(dy, y, p, z)
        dy[1] = y[2]
        dy[2] = y[3]
        dy[3] = y[4]
        dy[4] = -k2 * y[3] # From integration of X'''' + k^2 X'' = 0 -> X''' + k^2 X' = C (Shear)
                           # But strictly, the 4th order form is safer for general BVP solvers.
                           # Eq: EI*X'''' + T*X'' = 0
    end

    # Correct ODE for 4th order buckling:
    # Landau Eq (21.1): EI*X'''' + T*X'' = 0
    # Rearranged: X'''' = -(T/EI) * X''
    function rod_buckle_ode!(dy, y, p, z)
        dy[1] = y[2]       # X'
        dy[2] = y[3]       # X''
        dy[3] = y[4]       # X'''
        dy[4] = -k2 * y[3] # X'''' = -k^2 * X''
    end

    # CLAMPED BCs: Position and Slope are zero at both ends
    function bc_clamped!(residual, y, p, z)
        residual[1] = y[1][1]   # X(0) = 0
        residual[2] = y[1][2]   # X'(0) = 0 (No rotation)
        residual[3] = y[end][1] # X(L) = 0
        residual[4] = y[end][2] # X'(L) = 0 (No rotation)
    end

    zspan = (0.0, L)
    
    # Initial guess is tricky for eigenvalue problems. 
    # We need a slight asymmetry to trigger the mode or a determinant solver.
    # Here we solve for the shape given the KNOWN analytical T_crit for demonstration.
    # T_crit = 4 * pi^2 * EI / L^2
    
    prob = BVProblem(rod_buckle_ode!, bc_clamped!, [0.0, 0.0, 1.0, 0.0], zspan)
    sol = solve(prob, Shooting(Vern7()))
    return sol
end

# Parameters
L_val = 2.0
EI_val = 1000.0

# Critical Load for Clamped-Clamped is 4x that of Hinged
# T_hinged = pi^2 * EI / L^2
# T_clamped = 4 * pi^2 * EI / L^2
T_clamped = 4 * (pi^2 * EI_val) / L_val^2

sol_clamp = buckling_clamped(L_val, T_clamped, EI_val)

plot(sol_clamp, vars=(0, 1), 
     title="Buckling Mode (Clamped Ends)", 
     label="T_cr = 4π²EI/L²", lw=3, color=:blue,
     xlabel="z", ylabel="Deflection X")