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
    @views res[1:3] .= y[4:6] .- yd[1:3]
    @views res[4:6] .= yd[4:6] .- G_EARTH 
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
function plot(res::Result)
    plt.ax1 = plt.subplot(111) 
    plt.ax1.set_xlabel("time [s]")
    plt.plot(res.time, res.pos_z, color="green")
    plt.ax1.set_ylabel("pos_z [m]")  
    plt.ax1.grid(true) 
    plt.ax2 = plt.twinx()  
    plt.ax2.set_ylabel("vel_z [m/s]")   
    plt.plot(res.time, res.vel_z, color="red")
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
    s = nothing
    prob = DAEProblem(res!, yd0, y0, tspan, s, differential_vars=differential_vars)
    prob, solver
end
function solve!(res, prob, solver)
    abstol  = 0.0006 # max error in m/s and m
    for (i,t) in pairs(dt:dt:t_final)
        tspan2 = (t-dt, t)
        prob2   = remake(prob; tspan=tspan2)
        sol     = solve(prob2, solver; abstol, reltol=0.001)
        res.time[i]  = t
        res.pos_z[i] = sol.u[end][3]
        res.vel_z[i] = sol.u[end][6]
    end
    nothing
end
prob::DAEProblem, solver::IDA=init()
res = Result(2*t_final)
sol=solve!(res, prob, solver)

# @time time_, pos_z, vel_z = solve(integrator, dt, t_final)
# @time time_, pos_z, vel_z = solve(integrator, dt, 2*t_final)
nothing
