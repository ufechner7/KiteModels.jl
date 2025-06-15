## Testing the develop branch before merging it
On  AMD Ryzen 7 7840U.

### Pass criteria for the performance:
First run:
- simplifying the system      < 45s

Second run:
- system initialized          < 15s
- total time without plotting < 35s

## Current results on main, 15.06.2025:
```
julia> include("examples/ram_air_kite.jl")
[ Info: Loading packages 
Time elapsed: 2.769526784 s
[ Info: Creating wing, aero, vsm_solver, system_structure and symbolic_awe_model:
Time elapsed: 4.461972064 s
[ Info: Initialized integrator in 12.45877899 seconds
[ Info: System initialized at:
Time elapsed: 19.64399933 s
[ Info: Total time without plotting:
Time elapsed: 37.744739449 s
┌ Info: Performance:
│   times_realtime = 8.560932325389938
└   integrator_times_realtime = 57.896691251570246
```

### Develop branch
In Bash:
```
juliaup default 1.11
git checkout develop
cd bin
./update_default_manifest
./update_xz_file
cd ..
git commit -m "Update default manifest" Manifest-v1.11.toml.default
git push
jl
```
Now in Julia:
```
include("examples/ram_air_kite.jl")
```
Expected output:
```
julia> include("examples/ram_air_kite.jl")
[ Info: Loading packages 
Time elapsed: 2.863292943 s
[ Info: Creating wing, aero, vsm_solver, system_structure and s:
Time elapsed: 4.51880941 s
[ Info: Initialized integrator in 6.000656365 seconds
[ Info: System initialized at:
Time elapsed: 13.392478822 s
[ Info: Total time without plotting:
Time elapsed: 30.847251362 s
┌ Info: Performance:
│   times_realtime = 8.426564704784319
└   integrator_times_realtime = 62.47122497288717
```

### NLSolve test
```
include("mwes/mwe_26.jl")
include("mwes/mwe_26.jl")
include("mwes/mwe_26.jl")
```
Each time the force, that is printed must be > 9000 N.
