using DifferentialEquations
using Plots

function equations(dx,x,p,t)
    αₓ = 1.5
    βₓ = 5.0
    z_x = 0.4
    n_zx = 2.7
    x_z = 1.5
    n_xz = 2.7
    δ_z = 1
    S = p

 dx[1] = (αₓ + βₓ*S) / (1 + S + (x[2]/z_x)^n_zx) - x[1]
 dx[2] = 1 / (1 + (x[1]/x_z)^n_xz) - δ_z*x[2]
end

param = exp10.(range(-2 , stop = 2 , length = 5000))

X = zeros(size(param))
Z = zeros(size(param))

u0 = [0.2;1.5]
tspan = (0.0,20.0)

for i = 1:length(param)

    p = param[i]
    prob = ODEProblem(equations,u0,tspan,p)
    sol = solve(prob)

    X[i] = sol[1,end]
    Z[i] = sol[2,end]

end

plot(param,X,xaxis = ("S", :log),yaxis="X",label = "", linewidth=3,linecolor = :black)
#savefig("2c.pdf")
