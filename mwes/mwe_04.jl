using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using Sundials, ControlPlots

const G_EARTH  = [0.0, 0.0, -9.81] # gravitational acceleration
const dt = 0.05
const t_final = 10.0
  
# Example one: Falling mass
# State vector y   = mass.pos, mass.vel
# Derivative   yd  = mass.vel, mass.acc
# Residual     res = (y.vel - yd.vel), (yd.acc - G_EARTH)     
function res!(res, yd, y, p, time)
    res[1:3] .= y[4:6] .- yd[1:3]
    res[4:6] .= yd[4:6] .- G_EARTH 
end

struct Result
    time::Vector{Float64}
    pos_z::Vector{Float64}
    vel_z::Vector{Float64}
end
function Result(t_final)
    n=Int64(round(t_final/dt+1))
    Result(zeros(n), zeros(n), zeros(n))
end

function solve!(res, integrator, dt, t_final)
    for (i,t) in pairs(0:dt:t_final)
        res.time[i] = t
        Sundials.step!(integrator, dt, true)
        res.pos_z[i] = integrator.u[3]
        res.vel_z[i] = integrator.u[6]
    end
    nothing
end

function init()
    vel_0 = [0.0, 0.0, 50.0]    # Initial velocity
    pos_0 = [0.0, 0.0,  0.0]    # Initial position
    acc_0 = [0.0, 0.0, -9.81]   # Initial acceleration
    y0 = append!(pos_0, vel_0)  # Initial pos, vel
    yd0 = append!(vel_0, acc_0) # Initial vel, acc

    differential_vars = ones(Bool, length(y0))
    solver  = IDA(linear_solver=:GMRES, max_order = 4)
    tspan   = (0.0, t_final) 
    abstol  = 0.0006 # max error in m/s and m
    s = nothing

    prob    = DAEProblem(res!, yd0, y0, tspan, s, differential_vars=differential_vars)
    integrator = Sundials.init(prob, solver, abstol, reltol=0.001)
    return integrator
end

integrator=init()
res=Result(2*t_final)
@time solve!(res, integrator, dt, t_final)
@allocated solve!(res, integrator, dt, 2*t_final)

# 454864
