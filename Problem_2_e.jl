using DifferentialEquations
using Plots

function equations(dx,x,p,t)
    αₓ = 1.5 # 0.486
    α_y = 0.31 # 0.252
    βₓ = 7.17 # 5.31
    β_y = 14.9 # 6.2
    δ_y = 1.31 # 1.07
    δ_z = 1.41 # 2.0
    z_x = 1.75e-3 # 1.94e-3
    y_z = 0.13 # 0.097
    x_z = 0.19 # 0.038
    x_y = 2e-3  # 1.4e-3

    n_zx = 2 # 2.8

    n_xz = 2
    n_xy = 2.7
    n_yz = 2

    S = p

 dx[1] = (αₓ + βₓ*S) / (1 + S + (x[3]/z_x)^n_zx) - x[1]
 dx[2] = (α_y + β_y*S) / (1 + S + (x[1]/x_y)^n_xy) - δ_y*x[2]
 dx[3] = 1 / (1 + (x[1]/x_z)^n_xz +(x[2]/y_z)^n_yz) - δ_z*x[3]

end

# ----------- Below the bifurcation point ----------

p=80
# u0 = [0.5;10.0;1.0] works for p=10
u0 = [0.5;10.0;1.0]
tspan = (0.0,120.0)

prob_ss = SteadyStateProblem(equations,u0,p)
sol_ss = solve(prob_ss)
ss = Array(sol_ss)

#plot(sol)
#plot(sol,linewidth=3,xaxis="Time",yaxis="Concentration",label=["X" "Y" "Z"])
p=100
u0 = [ss[1];ss[2];ss[3]]
prob_1 = ODEProblem(equations,u0,tspan,p)
sol_1 = solve(prob_1)
plot(sol_1[3,:],linewidth=1,xaxis="Time",yaxis="Z",label="Cell 1")

p=100
u0 = 1.25*[ss[1];ss[2];ss[3]]
prob_2 = ODEProblem(equations,u0,tspan,p)
sol_2 = solve(prob_2)
plot!(sol_2[3,:],linewidth=1,xaxis="Time",yaxis="Z",label="Cell 2")

p=100
u0 = 0.75*[ss[1];ss[2];ss[3]]
prob_3 = ODEProblem(equations,u0,tspan,p)
sol_3 = solve(prob_3)
plot!(sol_3[3,:],linewidth=1,xaxis="Time",yaxis="Z",label="Cell 3")
#savefig("2e_hopf.png")

# --------- Above saddle node bifurcation point ------

p=5000
# u0 = [0.5;10.0;1.0] works for p=10
u0 = [10.0;10.0;1.0]
tspan = (0.0,120.0)

prob_ss = SteadyStateProblem(equations,u0,p)
sol_ss = solve(prob_ss)
ss = Array(sol_ss)

#plot(sol)
#plot(sol,linewidth=3,xaxis="Time",yaxis="Concentration",label=["X" "Y" "Z"])
p=100
u0 = [ss[1];ss[2];ss[3]]
prob_1 = ODEProblem(equations,u0,tspan,p)
sol_1 = solve(prob_1)
plot(sol_1[3,:],linewidth=1,xaxis="Time",yaxis="Z",label="Cell 1")
#plot(sol_1,linewidth=3,xaxis="Time",yaxis="Z",label="Cell 1")

p=100
u0 = 1.25*[ss[1];ss[2];ss[3]]
prob_2 = ODEProblem(equations,u0,tspan,p)
sol_2 = solve(prob_2)
plot!(sol_2[3,:],linewidth=1,xaxis="Time",yaxis="Z",label="Cell 2")
#plot(sol_1,linewidth=3,xaxis="Time",yaxis="Z",label="Cell 2")

p=100
u0 = 0.75*[ss[1];ss[2];ss[3]]
prob_3 = ODEProblem(equations,u0,tspan,p)
sol_3 = solve(prob_3)
plot!(sol_3[3,:],linewidth=1,xaxis="Time",yaxis="Z",label="Cell 3")
#plot(sol_1,linewidth=3,xaxis="Time",yaxis="Z",label="Cell 2")
#savefig("2e_saddle.png")
