using DifferentialEquations
using Plots

function equations(dx,x,p,t)
    αₓ = 0.039
    α_y = 0.0043
    βₓ = 6.1
    β_y = 5.7
    δ_y = 1.05
    δ_z = 1.04
    z_x = 1.3e-5
    y_z = 11e-3
    x_z = 0.12
    x_y = 7.9e-4
    n_zx = 2.32

    n_xz = 2
    n_xy = 2
    n_yz = 2

    S = p

 dx[1] = (αₓ + βₓ*S) / (1 + S + (x[3]/z_x)^n_zx) - x[1]
 dx[2] = (α_y + β_y*S) / (1 + S + (x[1]/x_y)^n_xy) - δ_y*x[2]
 dx[3] = 1 / (1 + (x[1]/x_z)^n_xz +(x[2]/y_z)^n_yz) - δ_z*x[3]

end

u0 = [0.0;0.0;0.0]
tspan = (0.0,50.0)

#plot(sol)
#plot(sol,linewidth=3,xaxis="Time",yaxis="Concentration",label=["X" "Y" "Z"])
p=0.02
prob = ODEProblem(equations,u0,tspan,p)
sol = solve(prob)
plot(sol[1,:],linewidth=3,xaxis="Time",yaxis="X",label="S = 0.02")
savefig("2d_0.02.png")

p=10.0
prob = ODEProblem(equations,u0,tspan,p)
sol = solve(prob)
plot(sol[1,:],linewidth=3,xaxis="Time",yaxis="X",label="S = 10")
savefig("2d_10.png")

p=10e5
prob = ODEProblem(equations,u0,tspan,p)
sol = solve(prob)
plot(sol[1,:],linewidth=3,xaxis="Time",yaxis="X",label="S = 10e5")
savefig("2d_10e5.png")
