using Revise, KiteModels, ModelingToolkit, ControlPlots

s = KPS4_3L(KCU(se("system_3l.yaml")))
update_settings()
s.set.segments = 2
s.set.aero_surfaces = 2
s.measure.winch_torque = [-0.2, -0.2, -70]
s.measure.tether_length = [s.set.l_tether, s.set.l_tether, s.set.l_tether]
s.measure.distance = s.measure.tether_length[3]
s.measure.elevation_left = deg2rad(85)
s.measure.elevation_right = deg2rad(85)
s.measure.azimuth_left = deg2rad(-5)
s.measure.azimuth_right = deg2rad(1)

prob, sol, ss = model!(s)

pos = prob[s.simple_sys.pos]
pos = [[pos[j, i] for j in 1:3] for i in 1:s.num_A]

l = s.set.l_tether+10
plot2d(pos, 0.0; zoom=true, front=false, xlim=(-l/2, l/2), ylim=(0, l))
# next_step!(s)
sol