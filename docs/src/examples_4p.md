```@meta
CurrentModule = KiteModels
```
# Examples for using the four point kite model

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
const kps = KPS4(kcu)
```
Then we call the function find_steady_state which uses a non-linear solver to find the solution for a given elevation angle, reel-out speed and wind speed. 
```julia
find_steady_state!(kps, prn=true)
```
Finding the steady state of the 4 point model is difficult and it only works when we artificially reduce the stiffness by a factor
of 0.035. In the function [`init_sim`](@ref) this factor is slowly increased to 1.0.

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
![Initial State](initial_state_4p.png)

## Print other model outputs
Print the vector of the positions of the particles:
```
julia> kps.pos
11-element StaticArrays.SVector{11, StaticArrays.MVector{3, Float64}} with indices SOneTo(11):
 [0.0, 0.0, 0.0]
 [27.956158502528176, 0.0, 61.20373231368065]
 [54.19504320716677, 0.0, 123.16303233303098]
 [78.43996882046399, 0.0, 185.92935677839677]
 [100.5072692264489, 0.0, 249.49394853600091]
 [120.24626566405448, 0.0, 313.81961163792636]
 [137.52679164398376, 0.0, 378.84882574376326]
 [138.4127462851352, 0.0, 383.77543513463166]
 [139.02795801127078, 0.0, 386.0060385324739]
 [138.8006044995987, 1.1208735809303805, 383.71585107239554]
 [138.8006044995987, -1.1208735809303805, 383.71585107239554]

```
Print the unstretched and stretched tether length and the height of the kite:
```julia
julia> unstretched_length(kps)
392.0

julia> tether_length(kps)
403.71695082721294

julia> calc_height(kps)
386.0060385324739
```
Because of the the stiffness_factor of 0.035 we have a longer tether-length then when using
the 1 point kite model. 

Print the force at the winch (groundstation, in Newton) and at each tether segment:
```julia
julia> winch_force(kps)
643.059283146209

julia> spring_forces(kps)
15-element Vector{Float64}:
 643.0081086309667
 642.9843324304267
 642.9695864010523
 642.9577694795955
 642.9475947472262
 642.938668170613
 165.89618777926103
 -25.722507476161205
  51.67720019071219
 254.61397516808276
 244.5988483585308
 244.5988483585308
 254.61397516808276
 -25.722507476161205
 166.98431814092817
```
Some of the forces are negative which means the segments are getting compressed. This is acceptable for
the kite itself (not for the tether).

Print the lift and drag forces of the kite (in Newton) and the lift over drag ratio:
```julia
julia> lift, drag = lift_drag(kps)
(616.7473148222452, 142.89285185868704)

julia> lift_over_drag(kps)
4.316152325325367
```
Print the wind speed vector at the kite:
```julia
julia> v_wind_kite(kps)
3-element StaticArrays.MVector{3, Float64} with indices SOneTo(3):
  13.310738776362681
  0.0
  0.0
```