# import Pkg; Pkg.add("DifferentialEquations")

using DifferentialEquations, Plots

#= 
-------------------------------------------------------------------------
 §20: Small Deflections (Numerical BVP)
 Rod with one end clamped and one end supported under uniform load q.
-------------------------------------------------------------------------
=#

function solve_rod_deflection(L, q, EI)
    # Define the 4th order ODE as a system of 1st order ODEs
    # y[1]=Y, y[2]=Y', y[3]=Y'', y[4]=Y'''
    function rod_ode!(dy, y, p, z)
        dy[1] = y[2]
        dy[2] = y[3]
        dy[3] = y[4]
        dy[4] = q / EI
    end

    # Boundary conditions: 
    # at z=0: Y=0, Y'=0
    # at z=L: Y=0, Y''=0
    function bc!(residual, y, p, z)
        residual[1] = y[1][1]       # Y(0) = 0
        residual[2] = y[1][2]       # Y'(0) = 0
        residual[3] = y[end][1]     # Y(L) = 0
        residual[4] = y[end][3]     # Y''(L) = 0
    end

    zspan = (0.0, L)
    prob = BVProblem(rod_ode!, bc!, [0.0, 0.0, 0.0, 0.0], zspan)
    sol = solve(prob, Shooting(Vern7()))
    return sol
end

# Parameters for a steel rod
L = 2.0; q = 100.0; E = 210e9; I = (pi * 0.05^4)/4
sol_deflect = solve_rod_deflection(L, q, E*I)

plot(sol_deflect, vars=(0, 1), title="Deflection of Clamped-Supported Rod",
     xlabel="z (Length)", ylabel="Y (Deflection)", label="Numerical Solution", lw=2)