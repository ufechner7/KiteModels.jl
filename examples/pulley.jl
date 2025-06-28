# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

using KiteModels, VortexStepMethod, ControlPlots

set = se("system_ram.yaml")
set.v_wind = 10.0
set.l_tether = 5.0
set.abs_tol = 1e-4
set.rel_tol = 1e-4
dynamics_type = DYNAMIC

points = Point[]
segments = Segment[]
pulleys = Pulley[]

push!(points, Point(1, [0.0, 0.0, 2.0], STATIC))
push!(points, Point(2, [2.0, 0.0, 2.0], STATIC))
push!(points, Point(3, [0.1, 0.0, 1.0], DYNAMIC))
push!(points, Point(4, [0.1, 0.0, 0.0], DYNAMIC; mass=0.1))

push!(segments, Segment(1, (3,1), BRIDLE))
push!(segments, Segment(2, (3,2), BRIDLE))
push!(segments, Segment(3, (3,4), BRIDLE))

push!(pulleys, Pulley(1, (1,2), DYNAMIC))

transforms = [Transform(1, -deg2rad(0.0), 0.0, 0.0; base_pos=[1.0, 0.0, 4.0], base_point_idx=1, rot_point_idx=2)]
sys_struct = KiteModels.SystemStructure("pulley", set; points, segments, pulleys, transforms)
plot(sys_struct, 0.0; zoom=false, l_tether=set.l_tether)

sam = SymbolicAWEModel(set, sys_struct)

init_sim!(sam; remake=false)
for i in 1:100
    plot(sam, i/set.sample_freq; zoom=false)
    next_step!(sam)
end
