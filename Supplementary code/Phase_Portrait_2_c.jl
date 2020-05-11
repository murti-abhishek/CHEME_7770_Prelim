include("PhasePortraitV2.jl")

# Function for a dual repression system without cooperativity
# x1: range of x1 values (i.e. A values)
# x2: range of x2 values (i.e. R values)
# We use `@.` to apply the calculations across all rows.
# Note that model parameters are specified within the function
# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)
function toggleMono(x1, x2)
    αₓ = 1.5
    βₓ = 5.0
    z_x = 0.4
    n_zx = 2.7
    x_z = 1.5
    n_xz = 2.7
    δ_z = 1
    S=2

    u = @. (αₓ + βₓ*S) / (1 + S + (x2/z_x)^n_zx) - x1
    v = @. 1 / (1 + (x1/x_z)^n_xz) - δ_z*x2

    return (u,v)
end

#Range of x1, x2 values
x1range = (0,10,15)          #Has the form (min, max, num points)
x2range = (0,2,15)          #Has the form (min, max, num points)
x₀ = ([2.0,2.0],[2.0, 0.5],[6.0, 0.2])  #initial state vectors; a common must be included after the first array
tspan=(0.0,30.0)             #time span

#Call the phaseplot functon to construct the phase portrait
phaseplot(toggleMono, x1range, x2range, xinit=x₀, t=tspan, clines=true,
        norm=true, scale=0.5)
