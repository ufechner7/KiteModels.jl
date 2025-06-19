# SPDX-FileCopyrightText: 2025 Bart van de Lint
#
# SPDX-License-Identifier: MPL-2.0

using KiteModels, VortexStepMethod, ControlPlots

set = se("system_ram.yaml")
set.l_tether = 5.0
dynamics_type = DYNAMIC

points = Point[]
segments = Segment[]
pulleys = Pulley[]

push!(points, Point(1, [0.0, 0.0, 2.0], STATIC))
push!(points, Point(2, [2.0, 0.0, 2.0], STATIC))
push!(points, Point(3, [1.0, 0.0, 1.0], DYNAMIC))
push!(points, Point(4, [1.1, 0.0, 0.0], DYNAMIC; mass=1.0))

push!(segments, Segment(1, (3,1), BRIDLE))
push!(segments, Segment(2, (3,2), BRIDLE))
push!(segments, Segment(3, (3,4), BRIDLE))

# push!(pulleys, Pulley(1, (1,2), DYNAMIC))

transforms = [Transform(1, 0.0, 0.0, 0.0; base_pos=[0.0, 0.0, 10.0], base_point_idx=1, rot_point_idx=2)]
sys_struct = KiteModels.SystemStructure("pulley", set; points, segments, pulleys, transforms)
plot(sys_struct, 0.0; zoom=false)

sam = SymbolicAWEModel(set, sys_struct)

init_sim!(sam; remake=true)
for i in 1:80
    plot(sam, i/set.sample_freq; zoom=false)
    next_step!(sam)
end
