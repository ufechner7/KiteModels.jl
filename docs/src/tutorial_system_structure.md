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

```julia
using KiteModels, VortexStepMethod, ControlPlots

set = se("system_ram.yaml")
set.segments = 20
dynamics_type = DYNAMIC
```
Then, we define vectors of the system structure types we are going to use. For this simple example we only need points and segments.
```julia
points = Point[]
segments = Segment[]

points = push!(points, Point(1, [0.0, 0.0, set.l_tether], STATIC; wing_idx=0))
```
The first point we add is a static point. There are four different [`DynamicsType`](@ref)s to choose from: `STATIC`, `QUASI_STATIC`, `DYNAMIC` and `WING`. `STATIC` just means that the point doesn't move. `DYNAMIC` is a point modeled with acceleration, while `QUASI_STATIC` constrains this acceleration to be zero at all times. A `WING` point is connected to a rigid wing body.

Now we can add `DYNAMIC` points and connect them to each other with segments. `BRIDLE` segments don't need to have a tether, because they have a constant unstretched length.
```julia
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
```julia
sys_struct = SystemStructure("tether"; points, segments)
plot(sys_struct, 0.0)
plt.gcf()
```

If the system looks good, we can easily model it, by first creating a [`SymbolicAWEModel`](@ref), initializing it and stepping through time.
```julia
model = SymbolicAWEModel(set, sys_struct)

init_sim!(model; remake=false)
for i in 1:100
    plot(model, i/set.sample_freq)
    next_step!(model; dt=1/set.sample_freq)
end
plt.gcf()
```
