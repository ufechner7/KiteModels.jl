```@meta
CurrentModule = KiteModels
```
# Custom SystemStructure and SymbolicAWESystem

A custom `SystemStructure` can be used to create models of kite power systems of almost any configuration.
- custom amount of tethers
- custom bridle configurations
- quasi-static or dynamic point masses
- different amounts of stiffness, damping and diameter on different tether segments

## Creating a simple tether

We start by loading the necessary packages and defining settings and parameters.

```@example 1
using ControlPlots
using KiteModels, VortexStepMethod

set = se("system_ram.yaml")
set.segments = 20
dynamics_type = DYNAMIC
```
Then, we define vectors of the system structure types we are going to use. For this simple example we only need points, that will be connected to eachother by segments.
```@example 1
points = Point[]
segments = Segment[]

points = push!(points, Point(1, [0.0, 0.0, set.l_tether], STATIC; wing_idx=0))
```
The first point we add is a static point. There are four different [`DynamicsType`](@ref)s to choose from: [`STATIC`](@ref), [`QUASI_STATIC`](@ref), [`DYNAMIC`](@ref) and [`WING`](@ref). `STATIC` just means that the point doesn't move. `DYNAMIC` is a point modeled with acceleration, while `QUASI_STATIC` constrains this acceleration to be zero at all times. A `WING` point is connected to a rigid wing body.

Now we can add `DYNAMIC` points and connect them to eachother with segments. `BRIDLE` segments don't need to have a tether, because they have a constant unstretched length.
```@example 1
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
```
From these arrays of points and segments we can create a [`SystemStructure`](@ref), which can be plotted in 2d to quickly investigate if the model is correct.
```@example 1
system_structure = SystemStructure("tether"; points, segments)
plot(system_structure, 0.0)
```

If the system looks good, we can easily model it, by first creating a [`SymbolicAWEModel`](@ref), initializing it and stepping through time.
```@example 1
model = SymbolicAWEModel(set, BodyAerodynamics[], VortexStepMethod.Solver[], system_structure)

init_sim!(model; remake=false)
plot(model, 0.0)
for i in 1:100
    next_step!(model; dt=1/set.sample_freq)
    plot(model, i/set.sample_freq)
end
```

## Creating a simple kite with one tether


