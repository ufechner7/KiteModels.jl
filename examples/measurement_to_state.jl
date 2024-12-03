using Revise, KiteModels, ModelingToolkit, ControlPlots, LinearAlgebra

function angle_between_vectors(v1, v2)
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2))
    θ = acos(clamp(cos_theta, -1.0, 1.0))  # Clamp to handle numerical errors
    return rad2deg(θ)  # Angle in radians
end

set = se("system_3l.yaml")
set.segments = 2
set.aero_surfaces = 2
s = KPS4_3L(KCU(set))
s.measure.winch_torque = [-0.2, -0.2, -70]
s.measure.tether_length = [47.37, 47.59, 47.639]
s.measure.distance = s.set.l_tether
s.measure.elevation_left = deg2rad(80)
s.measure.elevation_right = deg2rad(80)
s.measure.azimuth_left = deg2rad(-1)
s.measure.azimuth_right = deg2rad(1)
s.measure.distance_acc = s.measure.tether_acc[3]

prob, sol, ss, u0map = model!(s; real=false)

pos = sol[ss.pos]
pos = [[pos[j, i] for j in 1:3] for i in 1:s.num_A]

l = s.set.l_tether+10
# plot2d(pos, 0.0; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))

# next_step!(s)
@show angle_between_vectors(pos[5], pos[5+3]-pos[5])
@show sol.retcode
# for (u, e) in zip(sol.resid, equations(prob.f.sys))
#     println(u, "\t", e)
# end

t = 0
try
    while t < 20
        global t
        plot2d(s.pos, t; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))
        t = next_step!(s; set_values=[-0.2, -0.2, -70])
    end
catch AssertionError
    @show t
end
nothing