using Revise, KiteModels, OrdinaryDiffEqCore, OrdinaryDiffEqBDF, OrdinaryDiffEqSDIRK, LinearAlgebra, Statistics
using Base: summarysize

using Pkg
if ! ("ControlPlots" ∈ keys(Pkg.project().dependencies))
    using TestEnv; TestEnv.activate()
end
using ControlPlots

# TODO: sometimes very bad sim

set = deepcopy(load_settings("system_3l.yaml"))
# set.elevation = 71
# dt = 1.0
total_time = 5.0

steps = Int(round(total_time / dt))
logger = Logger(3*set.segments + 6, steps)

if !@isdefined s; s = KPS4_3L(KCU(set)); end
s.set = update_settings()
s.set.abs_tol = 0.0006
s.set.rel_tol = 0.001
s.set.l_tether = 50.0
# s.set.damping = 473
s.set.elevation = 60
init_set_values = [-0.1, -0.1, -200.0]
KiteModels.init_sim!(s; prn=true, torque_control=true, init_set_values, ϵ=10.0, flap_damping=0.1)
sys_state = KiteModels.SysState(s)


"""
Hypothesis: turn_rate_y = k * flap_diff
"""
diffs = zeros(1)
turn_rates = zeros(1)
steering = copy(init_set_values)
for dt in 0.5:0.1:2.0
    @show dt
    for i in 1:10
        steering[2] = -10 + i
        KiteModels.init_sim!(s; prn=false, torque_control=true, init_set_values, ϵ=0.1, flap_damping=0.1)
        KiteModels.next_step!(s; set_values=steering, dt=dt)
        s.integrator[s.simple_sys.heading_y]
        append!(diffs, s.get_flap_angle(s.integrator)[2]-s.get_flap_angle(s.integrator)[1])
        append!(turn_rates, s.integrator[s.simple_sys.heading_y]/dt)
        # l = s.set.l_tether+10
        # plot2d(s.pos, i-1; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))
    end
    plt.scatter(diffs, turn_rates)
    plt.xlabel("diffs")
    plt.ylabel("turn_rates")
    plt.show()
end
