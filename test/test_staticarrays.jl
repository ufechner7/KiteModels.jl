# unittest for StaticArrays
# fails for v1.5.13 and v1.5.14
using StaticArrays, LinearAlgebra, BenchmarkTools, Test

const KVec3    = MVector{3, Float64}

Base.@kwdef mutable struct KPS4{S, T, P}
    v_apparent::T =       zeros(S, 3)
end

const kps4 = KPS4{Float64, KVec3, 6+4+1}()

@inline function calc_particle_forces!(s, pos1, pos2)
    segment = pos1 - pos2
    norm1 = norm(segment)
    unit_vector = segment / norm1

    v_app_perp = s.v_apparent - s.v_apparent ⋅ unit_vector * unit_vector
    half_drag_force = norm(v_app_perp)
    nothing
end

# benchmark calc_particle_forces!
t = @benchmark calc_particle_forces!(kps4, pos1, pos2) setup=(pos1 = KVec3(1.0, 2.0, 3.0);  
                                        pos2 = KVec3(2.0, 3.0, 4.0))
@test t.memory == 0