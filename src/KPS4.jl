# implementation of a four point kite model
# to be included from KiteModels.jl

"""
    mutable struct KPS4{S, T, P} <: AbstractKiteModel

State of the kite power system, using a 4 point kite model. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. MVector{3, SimFloat}
- P: number of points of the system, segments+1

Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS3{S, T, P} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings = se()
    "Reference to the KCU struct (Kite Control Unit, type from the module KitePodSimulor"
    kcu::KCU = KCU()
    "Function for calculation the lift coefficent, using a spline based on the provided value pairs."
    calc_cl = Spline1D(se().alpha_cl, se().cl_list)
    "Function for calculation the drag coefficent, using a spline based on the provided value pairs."
    calc_cd = Spline1D(se().alpha_cd, se().cd_list)   
    "wind vector at the height of the kite" 
    v_wind::T =           zeros(S, 3)
    "wind vector at reference height" 
    v_wind_gnd::T =       zeros(S, 3)
    "wind vector used for the calculation of the tether drag"
    v_wind_tether::T =    zeros(S, 3)
    "apparent wind vector at the kite"
    v_apparent::T =       zeros(S, 3)
    v_app_perp::T =       zeros(S, 3)
    drag_force::T =       zeros(S, 3)
    lift_force::T =       zeros(S, 3)
    steering_force::T =   zeros(S, 3)
    last_force::T =       zeros(S, 3)
    spring_force::T =     zeros(S, 3)
    total_forces::T =     zeros(S, 3)
    force::T =            zeros(S, 3)
    unit_vector::T =      zeros(S, 3)
    av_vel::T =           zeros(S, 3)
    kite_y::T =           zeros(S, 3)
    segment::T =          zeros(S, 3)
    last_tether_drag::T = zeros(S, 3)
    acc::T =              zeros(S, 3)     
    vec_z::T =            zeros(S, 3)
    pos_kite::T =         zeros(S, 3)
    v_kite::T =           zeros(S, 3)        
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "area of one tether segment"
    seg_area::S =         zero(S) 
    bridle_area::S =      zero(S)
    "spring constant, depending on the length of the tether segment"
    c_spring::S =         zero(S)
    length::S =           0.0
    "damping factor, depending on the length of the tether segment"
    damping::S =          zero(S)
    area::S =             zero(S)
    last_v_app_norm_tether::S = zero(S)
    "lift coefficient of the kite, depending on the angle of attack"
    param_cl::S =         0.2
    "drag coefficient of the kite, depending on the angle of attack"
    param_cd::S =         1.0
    v_app_norm::S =       zero(S)
    cor_steering::S =     zero(S)
    "azimuth angle in radian; inital value is zero"
    psi::S =              zero(S)
    "elevation angle in radian; initial value about 70 degrees"
    beta::S =             deg2rad(se().elevation)
    last_alpha::S =        0.1
    alpha_depower::S =     0.0
    "relative start time of the current time interval"
    t_0::S =               0.0
    v_reel_out::S =        0.0
    last_v_reel_out::S =   0.0
    l_tether::S =          0.0
    rho::S =               0.0
    depower::S =           0.0
    steering::S =          0.0
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = ones(P)
end
