```@meta
CurrentModule = KiteModels
```
# Examples for using the one point kite model

## Create a test project
```bash
mkdir test
cd test
julia --project="."
```
and add the KiteUtils to the project:
```julia
]
add KiteUtils
<BACKSPACE>
```
Then, copy the default configuration files to your new project:
```julia
using KiteUtils
copy_settings()
```
Finally, add the KitePodModels and the KiteModels
```julia
]
add KitePodModels
add KiteModels
<BACKSPACE>
```

## Plotting the initial state
First an instance of the model of the kite control unit (KCU) is created which is needed by the Kite Power System model KPS3. Then we create a kps instance, passing the kcu model as parameter. We need to declare these variables as const to achieve a decent performance.
```julia
using KiteModels, KitePodModels, KiteUtils
const kcu = KCU(se())
const kps = KPS3(kcu)
```
Then we call the function find_steady_state which uses a non-linear solver to find the solution for a given elevation angle, reel-out speed and wind speed. 
```julia
find_steady_state!(kps, prn=true)
```
To plot the result in 2D we extract the vectors of the x and z coordinates of the tether particles with a for loop:
```julia
x = Float64[] 
z = Float64[]
for i in 1:length(kps.pos)
     push!(x, kps.pos[i][1])
     push!(z, kps.pos[i][3])
end
```
And finally we plot the postion of the particles in the x-z plane. When you type ```using Plots``` you will be ask if you want to install the Plots package. Just press \<ENTER\> and it gets installed.
```julia
using Plots
plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
plot!(x, z, seriestype = :scatter)
```
### Inital State
![Initial State](initial_state.png)

## Print other model outputs
Print the vector of the positions of the particles:
```
julia> kps.pos
7-element StaticArrays.SVector{7, StaticArrays.MVector{3, Float64}} with indices SOneTo(7):
 [0.0, 0.0, 0.0]
 [26.95751778658999, 0.0, 59.59749511924355]
 [51.97088814144287, 0.0, 120.03746888266994]
 [75.01423773175357, 0.0, 181.25637381120865]
 [96.06809940556136, 0.0, 243.18841293054678]
 [115.11959241520753, 0.0, 305.7661763854397]
 [132.79571663189674, 0.0, 368.74701279158705]

```
Print the unstretched and stretched tether length and the height of the kite:
```julia
julia> unstretched_length(kps)
392.0

julia> tether_length(kps)
392.4751313610764

julia> calc_height(kps)
368.74701279158705
``` 
Print the force at the winch (groundstation, in Newton) and at each tether segment:
```julia
julia> winch_force(kps)
728.5569144505084

julia> spring_forces(kps)
6-element Vector{Float64}:
 728.4835079763607
 734.9505623866943
 741.505320143339
 748.1408238767988
 754.8499002675924
 761.6993164647175
```
The force increases when going upwards because the kite not only experiances the winch force, but in addition the weight of the tether.

Print the lift and drag forces of the kite (in Newton) and the lift over drag ratio:
```julia
julia> lift, drag = lift_drag(kps)
(888.5715658243445, 188.25229350390242)

julia> lift_over_drag(kps)
4.720110173881757
```
Print the wind speed vector at the kite:
```julia
julia> v_wind_kite(kps)
3-element StaticArrays.MVector{3, Float64} with indices SOneTo(3):
  13.308227860928211
  0.0
  0.0
```
