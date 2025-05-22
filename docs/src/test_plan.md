## Testing the develop branch before merging it
On  AMD Ryzen 7 7840U.

### Pass criteria for the performance:
First run:
- simplifying the system      < 45s

Second run:
- system initialized          < 15s
- total time without plotting < 35s

### Develop branch
In Bash:
```
juliaup default 1.11
git checkout develop
cd bin
./update_default_manifest
cd ..
git commit -m "Update default manifest" Manifest-v1.11.toml.default
git push
rm data/prob_1.11*.bin
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
Time elapsed: 2.843903577 s
[ Info: Creating wing, aero, vsm_solver, point_system and s:
Time elapsed: 11.780494151 s
[ Info: Creating ODESystem
  4.189182 seconds (7.83 M allocations: 200.190 MiB, 0.69% gc time, 24.40% compilation time: 6% of which was recompilation)
[ Info: Simplifying the system
 38.041611 seconds (304.71 M allocations: 10.010 GiB, 2.12% gc time, 25.50% compilation time: 18% of which was recompilation)
[ Info: Creating ODEProblem
 84.739970 seconds (640.11 M allocations: 20.794 GiB, 1.83% gc time, 22.66% compilation time: 12% of which was recompilation)
[ Info: Initialized integrator in 22.253394034 seconds
[ Info: System initialized at:
Time elapsed: 193.93390819 s
[ Info: Total time without plotting:
Time elapsed: 211.095669066 s
┌ Info: Performance:
│   times_realtime = 8.017735204615859
└   integrator_times_realtime = 60.411279264296
```

### Main branch
```
julia> include("examples/ram_air_kite.jl")
[ Info: Loading packages 
Time elapsed: 2.9059458 s
[ Info: Creating wing, aero, vsm_solver, point_system and s:
Time elapsed: 11.790243385 s
┌ Warning: Failure to deserialize /home/ufechner/repos/KiteModels.jl/data/prob_dynamic_1.11_3_seg.bin !
└ @ KiteModels ~/repos/KiteModels.jl/src/ram_air_kite.jl:519
[ Info: Rebuilding the system. This can take some minutes...
[ Info: Creating ODESystem
  4.775993 seconds (8.31 M allocations: 201.174 MiB, 0.64% gc time, 20.63% compilation time: 6% of which was recompilation)
[ Info: Simplifying the system
 43.592286 seconds (329.71 M allocations: 10.926 GiB, 2.06% gc time, 25.35% compilation time: 18% of which was recompilation)
[ Info: Creating ODEProblem
 95.744396 seconds (682.85 M allocations: 23.818 GiB, 2.08% gc time, 17.67% compilation time: 14% of which was recompilation)
[ Info: Initialized integrator in 23.951264242 seconds
[ Info: System initialized at:
Time elapsed: 212.672152747 s
[ Info: Total time without plotting:
Time elapsed: 232.546116606 s
┌ Info: Performance:
│   times_realtime = 4.415094072484458
└   integrator_times_realtime = 14.752422253318498
```

### Update precompiled code
```
cd bin
./update_default_manifest
cd ..
cd data
cp prob_dynamic_1.11_3_seg.bin prob_dynamic_1.11_3_seg.bin.default
cd ..
echo " " >> src/precompile.jl 
```
## Test again
julia> include("examples/ram_air_kite.jl")
[ Info: Loading packages 
Time elapsed: 2.807525651 s
[ Info: Creating wing, aero, vsm_solver, point_system and s:
Time elapsed: 4.499499709 s
[ Info: Initialized integrator in 5.94955637 seconds
[ Info: System initialized at:
Time elapsed: 13.388212381 s
[ Info: Total time without plotting:
Time elapsed: 33.223129728 s
┌ Info: Performance:
│   times_realtime = 4.4603427235287265
└   integrator_times_realtime = 15.033714698785786