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
    n=Int64(round(t_final/dt)+1)
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

function init(res)
    vel_0 = [0.0, 0.0, 50.0]    # Initial velocity
    pos_0 = [0.0, 0.0,  0.0]    # Initial position
    acc_0 = [0.0, 0.0, -9.81]   # Initial acceleration
    y0 = append!(pos_0, vel_0)  # Initial pos, vel
    yd0 = append!(vel_0, acc_0) # Initial vel, acc
    res.pos_z[1] = y0[3]
    res.vel_z[1] = y0[6]
    differential_vars = ones(Bool, length(y0))
    solver  = IDA(linear_solver=:GMRES, max_order = 4)
    tspan   = (0.0, t_final) 
    s = nothing
    prob = DAEProblem(res!, yd0, y0, tspan, s, differential_vars=differential_vars)
    prob, solver
end
function solve!(res, prob, solver)
    local sol
    abstol  = 0.0006 # max error in m/s and m
    for (i,t) in pairs(dt:dt:t_final)
        tspan2 = (t-dt, t)
        prob2   = remake(prob; tspan=tspan2)
        if i > 1
            y0  =sol.u[end]
            yd0 =sol.du[end]
            prob2   = remake(prob2; u0=y0, du0=yd0)
        end
        sol = solve(prob2, solver; abstol, reltol=0.001)
        res.time[i+1]  = t
        res.pos_z[i+1] = sol.u[end][3]
        res.vel_z[i+1] = sol.u[end][6]
    end
    nothing
end
res = Result(t_final)
prob::DAEProblem, solver::IDA=init(res)
solve!(res, prob, solver)
res = Result(t_final)
prob::DAEProblem, solver::IDA=init(res)
bytes=@allocated solve!(res, prob, solver)
n=Int64(round(t_final/dt+1))
println("Allocated $(Int64(round(bytes/n))) bytes per iteration!")
# plot(res)

nothing
