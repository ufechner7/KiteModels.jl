using Test, BenchmarkTools, StaticArrays, LinearAlgebra, KiteUtils
using KiteModels, KitePodModels

if ! @isdefined SEGMENTS
    const SEGMENTS = se().segments
end
if ! @isdefined kcu
    const kcu = KCU(se())
    const kps = KPS3(kcu)
end

include("../src/consts.jl")

res1 = zeros(SVector{SEGMENTS+1, KiteModels.KVec3})
res2 = deepcopy(res1)
if ! @isdefined res3
    if USE_WINCH
        const res3 = vcat(reduce(vcat, vcat(res1, res2)), zeros(2))
    else
        const res3 = vcat(reduce(vcat, vcat(res1, res2)))
    end
end

msg=""
@testset verbose = true "KPS3 benchmarking....   " begin

function set_defaults()
    KiteModels.clear!(kps)
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
    KiteModels.clear!(kps)
    kps.set.l_tether = 392.0
    kps.set.elevation = 70.0
    kps.set.area = 10.0
    kps.set.rel_side_area = 50.0
    kps.set.v_wind = 9.1
    kps.set.mass = 6.2
    kps.set.c_s = 0.6
end

set_defaults()
function test_initial_condition(params::Vector)
    my_state = kps
    y0, yd0 = KiteModels.init(my_state, params)
    residual!(res3, yd0, y0, kps, 0.0)
    return norm(res3) # z component of force on all particles but the first
end

# Inputs:
# State vector state_y   = pos1, pos2, ..., posn, vel1, vel2, ..., veln
# Derivative   der_yd    = vel1, vel2, ..., veln, acc1, acc2, ..., accn
# Output:
# Residual     res = res1, res2 = pos1,  ..., vel1, ...
@testset "test_residual!       " begin
    res1 = zeros(SVector{SEGMENTS, KVec3})
    res2 = deepcopy(res1)
    res = reduce(vcat, vcat(res1, res2))
    if USE_WINCH
        res = vcat(res, zeros(2))
    end
    X = zeros(SimFloat, 2*kps.set.segments)
    y0, yd0 = KiteModels.init(kps, X; delta=1e-6)
    # println(y0)
    # println(yd0)
    p = kps
    t = 0.0
    clear!(kps)
    residual!(res, yd0, y0, p, t)
    res1 = res[1:3*SEGMENTS]
    res2 = res[3*SEGMENTS+1:end]
    @test res1 == zeros(3*(SEGMENTS))
    # TODO: add test for res2
    # println(res2)
end

t = @benchmark residual!(res, yd, y, p, t) setup = (res1 = zeros(SVector{SEGMENTS, KVec3}); res2 = deepcopy(res1); 
                                                               res = vcat(reduce(vcat, vcat(res1, res2)), zeros(2)); pos = deepcopy(res1);
                                                               pos[1] .= [1.0,2,3]; vel = deepcopy(res1); X = zeros(SimFloat, 2*kps.set.segments); 
                                                               (y0, yd0) = KiteModels.init(kps, X); yd=yd0; y=y0;
                                                               p = kps; t = 0.0)

@test t.memory <= 64
global msg = "Mean time residual! one point model: $(round(mean(t.times), digits=1)) ns"

end
println(msg)

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