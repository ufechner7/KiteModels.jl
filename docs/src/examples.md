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
```julia
using KiteModels
using KitePodSimulator
const kcu = KCU()
const kps = KPS3(kcu)
find_steady_state(kps, true)

x = Float64[] 
z = Float64[]
for i in 1:length(kps.pos)
     push!(x, kps.pos[i][1])
     push!(z, kps.pos[i][3])
end  

using Plots
plot(x,z, xlabel="x [m]", ylabel="z [m]", legend=false)
plot!(x, z, seriestype = :scatter)
```
### Inital State
![Initial State](initial_state.png)