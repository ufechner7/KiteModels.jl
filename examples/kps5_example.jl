using KiteModels
using Timers
using Pkg 
using DAEProblemLibrary
# if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
#     using TestEnv; TestEnv.activate()
# end
using ControlPlots
import OrdinaryDiffEqCore.init
import OrdinaryDiffEqCore.step!
import OrdinaryDiffEqCore.solve
import OrdinaryDiffEqCore
using Dierckx, Interpolations, Serialization, StaticArrays, LinearAlgebra, Statistics, Parameters, NLsolve,
      DocStringExtensions, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, NonlinearSolve, FiniteDiff, DifferentiationInterface
# using ModelingToolkit: Symbolics, @register_symbolic
# using ModelingToolkit, ControlPlots, Parameters
set = deepcopy(load_settings("system_kps5.yaml"))
kcu::KCU = KCU(set)
s = KPS5(kcu)

dt = 1/s.set.sample_freq
time_range = 0:dt:set.sim_time-dt
steps = length(time_range)
logger = Logger(KiteModels.points(s), steps)
#solver  = DImplicitEuler(autodiff=false)#@AutoFiniteDiff())
KiteModels.get_kite_points(s)
KiteModels.calc_initial_state(s)
KiteModels.init_sim!(s) #(integrate gen getter in initsim)
KiteModels.generate_getters!(s)
KiteModels.simulate(s, logger)
save_log(logger, "tmp")
lg = load_log("tmp")
function play(s, lg)
    dt = 1/s.set.sample_freq
    conn = KiteModels.getconnections(s)
    sl = lg.syslog
    total_segmentsvector = Vector{Int64}[]
    for conn_pair in conn
        push!(total_segmentsvector, Int64[conn_pair[1], conn_pair[2]])
    end
    # Add tether segments
    for i in 0:(s.set.segments-2)
        push!(total_segmentsvector, [6+i, 6+i+1])
    end
    # Add final connection from last tether point to bridle point
    push!(total_segmentsvector, [6+s.set.segments-1, 1])
    for step in 1:length(0:dt:s.set.sim_time)-1 #-s.set.dt
        # Get positions at this time step
        x = sl.X[step]
        y = sl.Y[step]
        z = sl.Z[step] 
        # Create points array for all points in the system
        pointsvector = Vector{Float64}[]
        for i in 1:KiteModels.points(s)                 #FIX THIS!
            push!(pointsvector, Float64[x[i], y[i], z[i]])
        end        
        # Calculate appropriate limits for the plot
        x_min, x_max = 0, 40
        z_min, z_max = 0, 60
        t = (dt) * (step-1)
        # Plot the kite system at this time step
        plot2d(pointsvector, total_segmentsvector, t;
               zoom = false,
               xlim = (x_min, x_max),
               ylim = (z_min, z_max)
        )
        # Add a small delay to control animation speed
        sleep(0.05)
    end
    nothing
end
play(s, lg)


# @with_kw mutable struct Settings_not_there_yet @deftype Float64
#     g_earth::Vector{Float64} = [0.0, 0.0, -9.81]           # gravitational acceleration [m/s²]
#     save::Bool = false                                     # save animation frames
#     # spring constants
#     springconstant_tether::Float64 = 614600.0              # TETHER unit spring constant [N]
#     springconstant_bridle::Float64 = 10000.0               # BRIDLE unit spring constant [N]
#     springconstant_kite::Float64 = 60000.0                 # KITE unit spring constant [N]
#     # damping
#     damping_tether::Float64 = 473                          # TETHER unit damping coefficient [Ns]
#     rel_damping_kite::Float64 = 6.0                        # KITE relative unit damping coefficient [-]
#     rel_damping_bridle:: Float64 = 6.0                     # BRIDLE relative unit damping coefficient [-]
# end
# @with_kw mutable struct KPS5
#     "Reference to the settings2 struct"
#     #set::Settings2 = Settings2()
#     sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
#     t_0::Float64 = 0.0
#     iter::Int64 = 0
#     prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
#     integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Nothing} = nothing
#     get_state::Function            = () -> nothing
# end

# # -----------------------------
# # Simulation Function
# # -----------------------------
# function simulate(s, logger)
#     dt = 1/s.set.sample_freq
#     tol = s.set.tol
#     tspan = (0.0, dt)
#     time_range = 0:dt:s.set.sim_time-dt
#     steps = length(time_range)
#     iter = 0
#     for i in 1:steps
#         next_step!(s; dt=s.set.dt)
#         u = s.get_state(s.integrator)
#         x = u[1][1, :]
#         y = u[1][2, :]
#         z = u[1][3, :]
#         iter += s.iter
#         sys_state = SysState{points(s)}()
#         sys_state.X .= x
#         sys_state.Y .= y
#         sys_state.Z .= z
#         println("iter: $iter", " steps: $steps")
#         log!(logger, sys_state)
#         println(x[end], ", ", sys_state.X[end])
#     end
#     println("iter: $iter", " steps: $steps")
#     return nothing
# end


function plot_front_view3(lg)
    display(plotxy(lg.y, lg.z;
    xlabel="pos_y [m]",
    ylabel="height [m]",
    fig="front_view"))
    nothing
end