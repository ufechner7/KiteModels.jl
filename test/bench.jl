using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils

using KiteModels, KitePodModels

const SEGMENTS = se().segments
if ! @isdefined kcu
    const kcu = KCU()
    const kps = KPS3(kcu)
end

function set_defaults()
    KiteModels.clear(kps)
    kps.set.l_tether = 150.0
    kps.set.elevation = 60.0
    kps.set.area = 20.0
    kps.set.rel_side_area = 50.0
    kps.set.v_wind = 8.0
    kps.set.mass = 11.4
    kps.set.damping =  2 * 473.0
    kps.set.alpha = 1.0/7
    kps.set.c_s = 0.6
end

function init_392()
    KiteModels.clear(kps)
    kps.set.l_tether = 392.0
    kps.set.elevation = 70.0
    kps.set.area = 10.0
    kps.set.rel_side_area = 50.0
    kps.set.v_wind = 9.1
    kps.set.mass = 6.2
    kps.set.c_s = 0.6
end

res1 = zeros(SVector{SEGMENTS+1, KiteModels.KVec3})
res2 = deepcopy(res1)
if ! @isdefined res3
    const res3 = reduce(vcat, vcat(res1, res2))
end

set_defaults()
function test_initial_condition(params::Vector)
    my_state = kps
    y0, yd0 = KiteModels.init(my_state, params)
    residual!(res3, yd0, y0, kps, 0.0)
    return norm(res3) # z component of force on all particles but the first
end

res = nothing
x= nothing
z= nothing
@testset "test_initial_residual" begin
    global res, x, z
    init_392()
    initial_x =  [-1.52505,  -3.67761,  -5.51761,  -6.08916,  -4.41371,  0.902124,  0.366393,  0.909132,  1.27537,  1.1538,  0.300657,  -1.51768]
    res=test_initial_condition(initial_x)

    my_state = kps
    kps.set.l_tether = 392.0
    kps.set.elevation = 70.0
    kps.set.area = 10.0
    kps.set.v_wind = 9.1
    kps.set.mass = 6.2
    KiteModels.clear(my_state)
    x = Float64[] 
    z = Float64[]
    for i in 1:length(my_state.pos)
        push!(x, my_state.pos[i][1])
        push!(z, my_state.pos[i][3])
    end  
end

@benchmark residual!(res, yd, y, p, t) setup = (res1 = zeros(SVector{SEGMENTS, KVec3}); res2 = deepcopy(res1); 
                                                               res = reduce(vcat, vcat(res1, res2)); pos = deepcopy(res1);
                                                               pos[1] .= [1.0,2,3]; vel = deepcopy(res1); y = reduce(vcat, vcat(pos, vel));
                                                               der_pos = deepcopy(res1); der_vel = deepcopy(res1); yd = reduce(vcat, vcat(der_pos, der_vel));
                                                               p = kps; t = 0.0)

# julia> include("test/bench.jl")
# [ Info: Precompiling KiteModels [b94af626-7959-4878-9336-2adc27959007]
# Test Summary:         |
# test_initial_residual | No tests
# BenchmarkTools.Trial: 10000 samples with 191 evaluations.
#  Range (min … max):  523.084 ns … 813.450 ns  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     524.743 ns               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   528.733 ns ±  20.955 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

#   █▄  ▁▁                                                        ▁
#   ███▄███▆▄▄▄▄▄▃▄▃▃▁▁▁▃▄▁▃▁▁▁▅▇▇▆▆▇▇▆▆▆▆▄▁▃▁▃▄▁▄▃▃▃▁▃▃▁▁▃▄▁▁▃▃▇ █
#   523 ns        Histogram: log(frequency) by time        666 ns <

#  Memory estimate: 0 bytes, allocs estimate: 0.

# julia> include("test/bench.jl")
# Test Summary:         |
# test_initial_residual | No tests
# BenchmarkTools.Trial: 10000 samples with 191 evaluations.
#  Range (min … max):  523.073 ns … 953.518 ns  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     524.859 ns               ┊ GC (median):    0.00%
#  Time  (mean ± σ):   530.646 ns ±  27.683 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

#   █▅  ▂▃                                                        ▁
#   ███▆███▆▆▆▅▁▅▄▄▅▄▅▃▄▃▃▃▃▁▄▇▇▇▆▇▇▇▆▇▆▅▄▃▅▄▃▁▃▁▄▃▄▁▄▄▁▁▃▁▁▃█▇▅▅ █
#   523 ns        Histogram: log(frequency) by time        676 ns <

#  Memory estimate: 0 bytes, allocs estimate: 0.