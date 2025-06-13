using KiteModels, VortexStepMethod, ControlPlots

set = se("system_ram.yaml")
set.segments = 20
dynamics_type = DYNAMIC

points = Point[]
segments = Segment[]

points = push!(points, Point(1, [0.0, 0.0, set.l_tether], STATIC; wing_idx=0))

segment_idxs = Int[]
for i in 1:set.segments
    global points, segments
    point_idx = i+1
    pos = [0.0, 0.0, set.l_tether] - set.l_tether / set.segments * [0.0, 0.0, i]
    push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
    segment_idx = i
    push!(segments, Segment(segment_idx, (point_idx-1, point_idx), BRIDLE))
    push!(segment_idxs, segment_idx)
end

system_structure = SystemStructure("tether"; points, segments)
plot(system_structure, 0.0)

model = SymbolicAWEModel(set, system_structure)

init_sim!(model; remake=false)
for i in 1:100
    plot(model, i/set.sample_freq)
    next_step!(model)
end