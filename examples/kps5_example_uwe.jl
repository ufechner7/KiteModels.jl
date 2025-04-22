using Timers
tic()
using KiteModels
using Pkg 
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots
toc()

include("plotting_kps5.jl")

set = deepcopy(load_settings("system_kps5.yaml"))
kcu::KCU = KCU(set)
s = KPS5(kcu)

dt = 1/s.set.sample_freq
time_range = 0:dt:set.sim_time-dt
steps = length(time_range)
logger = Logger(KiteModels.points(s), steps)
KiteModels.get_kite_points(s)
KiteModels.calc_initial_state(s)
init_sim!(s)
KiteModels.generate_getters!(s)
KiteModels.simulate(s, logger)
save_log(logger, "tmp")

function play(s, lg, front_view = false, side_view = true)
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
        
        # Use comparison operator == instead of assignment =
        if front_view == true
            for i in 1:KiteModels.points(s)
                # For front view: Y on horizontal axis, Z on vertical axis
                push!(pointsvector, Float64[y[i], z[i], x[i]])
            end
            
            # Calculate appropriate limits for front view
            y_min, y_max = -10, 10
            z_min, z_max = 0, 50
            t = (dt) * (step-1)
            
            # Plot the kite system at this time step (front view)
            plot2d(pointsvector, total_segmentsvector, t;
                   zoom = false,
                   xlim = (y_min, y_max),
                   ylim = (z_min, z_max)
            )
        elseif side_view == true
            for i in 1:KiteModels.points(s)
                # For side view: X on horizontal axis, Z on vertical axis
                push!(pointsvector, Float64[x[i], y[i], z[i]])
            end
            
            # Calculate appropriate limits for side view
            x_min, x_max = 0, 20
            z_min, z_max = 0, 20
            t = (dt) * (step-1)
            
            # Plot the kite system at this time step (side view)
            plot2d(pointsvector, total_segmentsvector, t;
                   zoom = false,
                   xlim = (x_min, x_max),
                   ylim = (z_min, z_max)
            )
        end
        
        # Add a small delay to control animation speed
        sleep(0.05)
    end
    nothing
end

lg = load_log("tmp")
play(s, lg)