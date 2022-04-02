```@meta
CurrentModule = KiteModels
```
# Examples

## Create a test project
```bash
mkdir test
cd test
julia --project
```
and add KiteModels to the project:
```julia
]activate .
add KiteUtils
add KitePodSimulator
add KiteModels
<BACKSPACE>
```
finally, copy the default configuration files to your new project:
```julia
using KiteUtils
copy_settings()
```

## Plotting the initial state
First an instance of the model of the kite control unit (KCU) is created which is needed by the Kite Power System model KPS3. Then we create a kps instance, passing the kcu model as parameter. We need to declare these variables as const to achieve a decent performance.
```julia
using KiteModels
using KitePodModels
const kcu = KCU()
const kps = KPS3(kcu)
```
Then we call the function find_steady_state which uses a non-linear solver to find the solution for a given elevation angle, reel-out speed and wind speed. 
```julia
find_steady_state(kps, true)
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
 [0.0, 1.0e-6, 0.0]
 [26.957523083220014, 1.0e-6, 59.597492705934215]
 [51.97089867461361, 1.0e-6, 120.03746428611659]
 [75.01425347614648, 1.0e-6, 181.25636723731242]
 [96.06812033613873, 1.0e-6, 243.18840457807227]
 [115.11961850809006, 1.0e-6, 305.76616644603797]
 [132.79574786242222, 1.0e-6, 368.747001395149]
```
Print the unstretched and and stretched tether length:
```julia
julia> unstretched_length(kps)
392.0

julia> tether_length(kps)
392.4751313610764
``` 
Print the force at the winch (groundstation, in Newton) and at each tether segment:
```julia
julia> winch_force(kps)
728.5567740002092

julia> spring_forces(kps)
6-element Vector{Float64}:
 728.4833579505422
 734.950422647022
 741.5051811137938
 748.1406855651342
 754.8497626815621
 761.6991795967015
```
The force increases when going upwards because the kite not only experiances the winch force, but in addition the weight of the tether.

Print the lift and drag forces of the kite (in Newton) and the lift over drag ratio:
```julia
julia> lift, drag = lift_drag(kps)
(888.5714473490408, 188.25226817881344)

julia> lift_over_drag(kps)
4.720110179522627
```
Print the wind speed vector at the kite:
```julia
julia> v_wind_kite(kps)
3-element StaticArrays.MVector{3, Float64} with indices SOneTo(3):
 13.308227837486344
  0.0
  0.0
```

