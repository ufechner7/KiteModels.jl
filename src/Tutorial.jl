using StaticArrays, Parameters, BenchmarkTools, LinearAlgebra

@with_kw mutable struct KPS{T}
    "wind vector used for the calculation of the tether drag"
    v_wind_tether::T =    zeros(3)
    "apparent wind vector"
    v_apparent::T =       zeros(3)
    "vector of the drag force of the last calculated tether segment"
    last_tether_drag::T = zeros(3)  
end

const kps=KPS{MVector{3, Float64}}()
const kps2 = KPS{Vector{Float64}}()

function calc_drag(s::KPS, v_segment, unit_vector, rho, v_app_perp, area)
    s.v_apparent .= s.v_wind_tether - v_segment
    v_app_norm = norm(s.v_apparent)
    v_app_perp .= s.v_apparent .- s.v_apparent ⋅ unit_vector .* unit_vector
    s.last_tether_drag .= -0.5 * rho * norm(v_app_perp) * area .* v_app_perp
    v_app_norm
end 

@benchmark calc_drag($kps2, v_segment, unit_vector, rho, v_app_perp, area) setup=(v_segment=zeros(3); unit_vector=[0,1.0,0]; rho=1.0; v_app_perp=[10.0,0,0]; area=0.3)
# 155ns no in-place (dot) operations, dynamic vectors
#  70ns in-place operations, dynamic vectors
#  20ns in-place operations, static vectors (MVector)
