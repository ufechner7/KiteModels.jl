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
using KitePodSimulator
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
And finally we plot the result. When you type ```using Plots``` you will be ask if you want to install the Plots package. Just press \<ENTER\> and it gets installed.
```julia
using Plots
plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
plot!(x, z, seriestype = :scatter)
```
### Inital State
![Initial State](initial_state.png)