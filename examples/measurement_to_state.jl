using Revise, KiteModels, ModelingToolkit, ControlPlots, LinearAlgebra

function angle_between_vectors(v1, v2)
    cos_theta = dot(v1, v2) / (norm(v1) * norm(v2))
    θ = acos(clamp(cos_theta, -1.0, 1.0))  # Clamp to handle numerical errors
    return rad2deg(θ)  # Angle in radians
end

s = KPS4_3L(KCU(se("system_3l.yaml")))
update_settings()
s.set.segments = 2
s.set.aero_surfaces = 2
s.measure.winch_torque = [-0.2, -0.2, -20]
s.measure.tether_length = [s.set.l_tether, s.set.l_tether, s.set.l_tether]
s.measure.distance = s.measure.tether_length[3]
s.measure.elevation_left = deg2rad(89)
s.measure.elevation_right = deg2rad(89)
s.measure.azimuth_left = deg2rad(-3)
s.measure.azimuth_right = deg2rad(1)
s.measure.distance_acc = s.measure.tether_acc[3]

prob, sol, ss = model!(s)

pos = prob[s.simple_sys.pos]
pos = [[pos[j, i] for j in 1:3] for i in 1:s.num_A]

l = s.set.l_tether+10
plot2d(pos, 0.0; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))

# next_step!(s)
@show angle_between_vectors(s.pos[5], s.pos[5+3]-s.pos[5])
@show sol.retcode
nothing