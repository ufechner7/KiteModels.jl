# plot the lift and drag coefficients as function of angle of attack

using Printf
using KiteModels, KitePodModels, KiteUtils, LinearAlgebra

set = deepcopy(load_settings("system_v9.yaml"))

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots
plt.close("all")

set.abs_tol=0.00006
set.rel_tol=0.000001
V_WIND = 14.5

# the following values can be changed to match your interest
dt = 0.05
set.solver="DFBDF" # IDA or DFBDF
STEPS = 500
PLOT = true
PRINT = true
STATISTIC = false
DEPOWER = 0.47:-0.005:0.355
# end of user parameter section #

bridle_length = KiteModels.bridle_length(set)
println("bridle_length: $bridle_length")
bridle_area = (set.d_line/2000) * bridle_length
println("bridle_area: $bridle_area")

function set_tether_diameter!(se, d; c_spring_4mm = 614600, damping_4mm = 473)
    set.d_tether = d
    set.c_spring = c_spring_4mm * (d/4.0)^2
    set.damping = damping_4mm * (d/4.0)^2
end

set_tether_diameter!(set, set.d_tether)

elev = set.elevation
i = 1
set.v_wind = V_WIND # 25
logger = Logger(set.segments + 5, STEPS)
# set.depower = 100*depower
# set.depower_gain = 5

kcu = KCU(set)
kps4 = KPS4(kcu)
integrator = KiteModels.init_sim!(kps4; delta=0.03, stiffness_factor=0.5, prn=STATISTIC)
lift, drag = lift_drag(kps4)
println("Lift: $lift, Drag: $drag, Iterations: $(kps4.iter)")
