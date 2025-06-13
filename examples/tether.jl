using KiteModels, VortexStepMethod

include("plotting.jl")
set = se("system_ram.yaml")

points = Point[]
segments = Segment[]
dynamics_type = DYNAMIC

# First, create a static ground point that doesn't move
points = push!(points, Point(1, zeros(3), STATIC; wing_idx=0))

# Now, create the tether points
segment_idxs = Int[]
for i in 2:10
    global points, segments
    point_idx = i
    pos = [0.0, 0.0, -(i-1)]
    push!(points, Point(point_idx, pos, dynamics_type; wing_idx=0))
    segment_idx = i-1
    push!(segments, Segment(segment_idx, (point_idx-1, point_idx), POWER))
    push!(segment_idxs, segment_idx)
end
tethers = [Tether(1, segment_idxs)]

system_structure = SystemStructure("tether"; points, segments, tethers)
model = SymbolicAWEModel(set, BodyAerodynamics[], VortexStepMethod.Solver[], system_structure)

init_sim!(model)
next_step!(model)
