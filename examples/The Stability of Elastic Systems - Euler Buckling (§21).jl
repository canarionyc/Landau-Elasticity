#= 
-------------------------------------------------------------------------
 §21 Prob 1: Critical Compression (Buckling)
 Finding the first buckling mode numerically.
-------------------------------------------------------------------------
=#

function buckling_mode(L, T_val, EI)
    # Normalized ODE: X'''' + (T/EI)X'' = 0
    k2 = T_val / EI
    function buck_ode!(dy, y, p, z)
        dy[1] = y[2]
        dy[2] = y[3]
        dy[3] = y[4]
        dy[4] = -k2 * y[3]
    end

    # Boundary conditions for hinged ends: X(0)=X''(0)=0, X(L)=X''(L)=0
    function buck_bc!(residual, y, p, z)
        residual[1] = y[1][1]
        residual[2] = y[1][3]
        residual[3] = y[end][1]
        residual[4] = y[end][3]
    end

    zspan = (0.0, L)
    # We provide a non-zero initial guess to find the non-trivial solution
    prob = BVProblem(buck_ode!, buck_bc!, [0.0, 1.0, 0.0, -k2], zspan)
    sol = solve(prob, Shooting(Vern7()))
    return sol
end

# Analytical critical force
T_crit = (pi^2 * E * I) / L^2
sol_buckle = buckling_mode(L, T_crit, E*I)

plot(sol_buckle, vars=(0, 1), title="First Buckling Mode (Euler Buckling)",
     xlabel="z", ylabel="Lateral Displacement X", label="Mode 1", color=:red)