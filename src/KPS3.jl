#= MIT License

Copyright (c) 2020, 2021, 2022 Uwe Fechner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. =#

#= Model of a kite-power system in implicit form: residual = f(y, yd)

This model implements a 3D mass-spring system with reel-out. It uses five tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled as additional mass at the end of the tether.
The spring constant and the damping decrease with the segment length. The aerodynamic kite forces are
calculated, depending on reel-out speed, depower and steering settings. 

One point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#

"""
    mutable struct KPS3{S, T, P} <: AbstractKiteModel

State of the kite power system. Parameters:
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
    "Iteration"
    iter:: Int64 = 0
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
    "vector, perpendicular to v_apparent; output of calc_drag"
    v_app_perp::T =       zeros(S, 3)
    "drag force of kite and bridle; output of calc_aero_forces"
    drag_force::T =       zeros(S, 3)
    "lift force of the kite; output of calc_aero_forces"
    lift_force::T =       zeros(S, 3)
    "steering force acting on the kite; output of calc_aero_forces"
    steering_force::T =   zeros(S, 3)
    last_force::T =       zeros(S, 3)
    "spring force of the current tether segment, output of calc_res"
    spring_force::T =     zeros(S, 3)
    total_forces::T =     zeros(S, 3)
    "sum of spring and drag forces acting on the current segment, output of calc_res"
    force::T =            zeros(S, 3)
    "unit vector in the direction of the current tether segment, output of calc_res"
    unit_vector::T =      zeros(S, 3)
    "average velocity of the current tether segment, output of calc_res"
    av_vel::T =           zeros(S, 3)
    "y-vector of the kite fixed referense frame, output of calc_aero_forces"
    kite_y::T =           zeros(S, 3)
    "vector representing one tether segment (p1-p2)"
    segment::T =          zeros(S, 3)
    "vector of the drag force of the last calculated tether segment"
    last_tether_drag::T = zeros(S, 3)  
    vec_z::T =            zeros(S, 3) 
    "part one of the residual, difference between pos' and vel, non-flat, mainly for unit testing"
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "part two of the residual, difference between vel' and acc, non-flat, mainly for unit testing"
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "vector of the positions of the particles"
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "area of one tether segment"
    seg_area::S =         zero(S) 
    bridle_area::S =      zero(S)
    "spring constant, depending on the length of the tether segment"
    c_spring::S =         zero(S)
    "unstressed segment length [m]"
    segment_length::S =           0.0
    "damping factor, depending on the length of the tether segment"
    damping::S =          zero(S)
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
    "reel out speed of the winch [m/s]"
    v_reel_out::S =        0.0
    "reel out speed during the last time step"
    last_v_reel_out::S =   0.0
    "unstretched tether length"
    l_tether::S =          0.0
    "air density at the height of the kite"
    rho::S =               0.0
    depower::S =           0.0
    steering::S =          0.0
    "factor for the tether stiffness, used to find the steady state with a low stiffness first"
    stiffness_factor::S =  1.0
    "pre-calculated constant for the wind profile law calcuation"
    log_href_over_z0::S =  log(se().h_ref / se().z0)
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = ones(P)
end

"""
    clear(s::KPS3)

Initialize the kite power model.
"""
function clear(s::KPS3)
    s.t_0 = 0.0                              # relative start time of the current time interval
    s.v_reel_out = 0.0
    s.last_v_reel_out = 0.0
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    height = sin(deg2rad(s.set.elevation)) * s.set.l_tether
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s, height)
    s.alpha_depower = 0.0
    s.l_tether = s.set.l_tether
    s.segment_length = s.l_tether / s.set.segments
    s.beta = deg2rad(s.set.elevation)
    s.rho = s.set.rho_0
    mass_per_meter = s.set.rho_tether * π * (s.set.d_tether/2000.0)^2
    mass_per_meter = 0.011
    s.initial_masses .= ones(s.set.segments+1) * mass_per_meter * s.set.l_tether / s.set.segments # Dyneema: 1.1 kg/ 100m
    s.c_spring = s.set.c_spring / s.segment_length
    s.damping  = s.set.damping / s.segment_length
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list) 
end

function KPS3(kcu::KCU)
    s = KPS3{SimFloat, KVec3, kcu.set.segments+1}()
    s.set = kcu.set
    s.kcu = kcu
    s.calc_cl = Spline1D(s.set.alpha_cl, s.set.cl_list)
    s.calc_cd = Spline1D(s.set.alpha_cd, s.set.cd_list)       
    clear(s)
    return s
end

"""
    calc_drag(s::KPS3, v_segment, unit_vector, rho, last_tether_drag, v_app_perp)

Calculate the drag of one tether segment, result stored in parameter `last_tether_drag`.
Return the norm of the apparent wind velocity.
"""
function calc_drag(s::KPS3, v_segment, unit_vector, rho, v_app_perp, area)
    s.v_apparent .= s.v_wind_tether - v_segment
    v_app_norm = norm(s.v_apparent)
    v_app_perp .= s.v_apparent .- dot(s.v_apparent, unit_vector) .* unit_vector
    s.last_tether_drag .= -0.5 * s.set.cd_tether * rho * norm(v_app_perp) * area .* v_app_perp
    v_app_norm
end 

#     pos_kite:     position of the kite
#     rho:          air density [kg/m^3]
#     paramCD:      drag coefficient (function of power settings)
#     paramCL:      lift coefficient (function of power settings)
#     rel_steering: value between -1.0 and +1.0
function calc_aero_forces(s::KPS3, pos_kite, v_kite, rho, rel_steering)
    s.v_apparent    .= s.v_wind - v_kite
    s.v_app_norm     = norm(s.v_apparent)
    s.drag_force    .= s.v_apparent ./ s.v_app_norm
    s.kite_y        .= normalize(cross(pos_kite, s.drag_force))
    K                = 0.5 * rho * s.v_app_norm^2 * s.set.area
    s.lift_force    .= K * s.param_cl .* normalize(cross(s.drag_force, s.kite_y))   
    # some additional drag is created while steering
    s.drag_force    .*= K * s.param_cd * BRIDLE_DRAG * (1.0 + 0.6 * abs(rel_steering)) 
    s.cor_steering    = s.set.c2_cor / s.v_app_norm * sin(s.psi) * cos(s.beta) # in paper named i_(s,c), Eq. 30
    s.steering_force .= -K * s.set.rel_side_area/100.0 * s.set.c_s * (rel_steering + s.cor_steering) .* s.kite_y
    s.last_force     .= -(s.lift_force + s.drag_force + s.steering_force) 
    nothing
end

"""
    calc_height(s::KPS3)

Determine the height of the kite particle above ground.
"""
function calc_height(s::KPS3)
    pos_kite = s.pos[end]
    pos_kite[3]
end

# Calculate the vector res1, that depends on the velocity and the acceleration.
# The drag force of each segment is distributed equaly on both particles.
function calc_res(s::KPS3, pos1, pos2, vel1, vel2, mass, veld, result, i)
    s.segment .= pos1 - pos2
    height = (pos1[3] + pos2[3]) * 0.5
    rho = calc_rho(s, height)               # calculate the air density
    rel_vel = vel1 - vel2                # calculate the relative velocity
    s.av_vel .= 0.5 * (vel1 + vel2)
    norm1 = norm(s.segment)
    s.unit_vector .= normalize(s.segment) # unit vector in the direction of the tether
    # # look at: http://en.wikipedia.org/wiki/Vector_projection
    # # calculate the relative velocity in the direction of the spring (=segment)
    spring_vel = dot(s.unit_vector, rel_vel)

    k2 = 0.05 * s.c_spring * s.stiffness_factor             # compression stiffness tether segments
    if norm1 - s.segment_length > 0.0
        s.spring_force .= (s.c_spring * s.stiffness_factor * (norm1 - s.segment_length) + s.damping * spring_vel) .* s.unit_vector
    else
        s.spring_force .= k2 * ((norm1 - s.segment_length) + (s.damping * spring_vel)) .* s.unit_vector
    end
    s.seg_area = norm1 * s.set.d_tether/1000.0
    s.last_v_app_norm_tether = calc_drag(s, s.av_vel, s.unit_vector, rho, s.v_app_perp, s.seg_area)
    s.force .= s.spring_force + 0.5 * s.last_tether_drag

    if i == s.set.segments+1 # add the drag of the bridle lines
        s.bridle_area =  s.set.l_bridle * s.set.d_line/1000.0
        s.last_v_app_norm_tether = calc_drag(s, s.av_vel, s.unit_vector, rho, s.v_app_perp, s.bridle_area)
        s.force .+= s.last_tether_drag  
    end
   
    s.total_forces .= s.force + s.last_force
    s.last_force .= 0.5 * s.last_tether_drag - s.spring_force
    acc = s.total_forces ./ mass # create the vector of the spring acceleration
    # result .= veld - (s.acc + SVector(0,0, -G_EARTH)) # Python code, wrong
    result .= veld - (SVector(0, 0, -G_EARTH) - acc)
    nothing
end

# Calculate the vector res1 using a vector expression, and calculate res2 using a loop
# that iterates over all tether segments. 
function loop(s::KPS3, pos, vel, posd, veld, res1, res2)
    s.masses               .= s.segment_length / (s.l_tether / s.set.segments) .* s.initial_masses
    s.masses[s.set.segments+1]   += (s.set.mass + s.set.kcu_mass)
    res1[1] .= pos[1]
    res2[1] .= vel[1]
    for i in 2:s.set.segments+1
        res1[i] .= vel[i] - posd[i]
    end
    for i in s.set.segments+1:-1:2
        calc_res(s, pos[i], pos[i-1], vel[i], vel[i-1], s.masses[i], veld[i],  res2[i], i)
    end
    nothing
end

"""
    residual!(res, yd, y::MVector{S, SimFloat}, s::KPS3, time) where S

    N-point tether model, one point kite at the top:
    Inputs:
    State vector y   = pos1, pos2, ..., posn, vel1, vel2, ..., veln
    Derivative   yd  = vel1, vel2, ..., veln, acc1, acc2, ..., accn
    Output:
    Residual     res = res1, res2 = pos1,  ..., vel1, ...

    Additional parameters:
    s: Struct with work variables, type KPS3
    S: The dimension of the state vector
The number of the point masses of the model N = S/6, the state of each point 
is represented by two 3 element vectors.
"""
function residual!(res, yd, y::MVector{S, SimFloat}, s::KPS3, time) where S
    # unpack the vectors y and yd
    part = reshape(SVector{S}(y),  Size(3, div(S,6), 2))
    partd = reshape(SVector{S}(yd),  Size(3, div(S,6), 2))
    pos1, vel1 = part[:,:,1], part[:,:,2]
    pos = SVector{div(S,6)+1}(if i==1 SVector(0.0,0,0) else SVector(pos1[:,i-1]) end for i in 1:div(S,6)+1)
    vel = SVector{div(S,6)+1}(if i==1 SVector(0.0,0,0) else SVector(vel1[:,i-1]) end for i in 1:div(S,6)+1)
    posd1, veld1 = partd[:,:,1], partd[:,:,2]
    posd = SVector{div(S,6)+1}(if i==1 SVector(0.0,0,0) else SVector(posd1[:,i-1]) end for i in 1:div(S,6)+1)
    veld = SVector{div(S,6)+1}(if i==1 SVector(0.0,0,0) else SVector(veld1[:,i-1]) end for i in 1:div(S,6)+1)

    # update parameters
    pos_kite = pos[div(S,6)+1]
    v_kite   = vel[div(S,6)+1]
    delta_t = time - s.t_0
    delta_v = s.v_reel_out - s.last_v_reel_out
    s.segment_length = (s.l_tether + s.last_v_reel_out * delta_t + 0.5 * delta_v * delta_t^2) / div(S,6)
    s.c_spring = s.set.c_spring / s.segment_length
    s.damping  = s.set.damping / s.segment_length

    # call core calculation routines
    vec_c = SVector{3, SimFloat}(pos[s.set.segments] - pos_kite)     # convert to SVector to avoid allocations
    v_app = SVector{3, SimFloat}(s.v_wind - v_kite)
    calc_set_cl_cd(s, vec_c, v_app)
    calc_aero_forces(s, pos_kite, v_kite, s.rho, s.steering) # force at the kite
    loop(s, pos, vel, posd, veld, s.res1, s.res2)
  
    # copy and flatten result
    for i in 2:div(S,6)+1
        for j in 1:3
           @inbounds res[3*(i-2)+j]              = s.res1[i][j]
           @inbounds res[3*(div(S,6))+3*(i-2)+j] = s.res2[i][j]
        end
    end
    if norm(res) < 1e5
        # println(norm(res))
        for i in 1:length(pos)
            @inbounds s.pos[i] .= pos[i]
        end
    end
    # @assert ! isnan(norm(res))
    s.iter += 1
    nothing
end

# Calculate the initial conditions y0 and yd0. Tether with the initial elevation angle
# se().elevation, particle zero fixed at origin.
# Parameters:
# x: vector of deviations of the tether particle positions from a straight line in x and z
# length(x) == 2*SEGMENTS
# Returns:
# res, a single vector consisting of the elements of y0 and yd0
function init_inner(s::KPS3, X=zeros(2 * s.set.segments); old=false, delta=0.0)
    pos = zeros(SVector{s.set.segments+1, KVec3})
    vel = zeros(SVector{s.set.segments+1, KVec3})
    acc = zeros(SVector{s.set.segments+1, KVec3})
    state_y0 = zeros(SVector{2*s.set.segments, KVec3})
    yd0 = zeros(SVector{2*s.set.segments, KVec3})

    DELTA = delta
    set_cl_cd(s, 10.0/180.0 * π)

    for i in 0:s.set.segments
        radius =  -i * s.set.l_tether / s.set.segments
        elevation = s.set.elevation
        sin_el, cos_el = sin(elevation / 180.0 * π), cos(elevation / 180.0 * π)
        if i == 0
            pos[i+1] .= SVec3(0.0, DELTA, 0.0)
        else
            pos[i+1] .= SVec3(-cos_el * radius+X[i], DELTA, -sin_el * radius+X[s.set.segments+i])
        end
        vel[i+1] .= SVec3(DELTA, DELTA, DELTA)
        acc[i+1] .= SVec3(DELTA, DELTA, DELTA)
    end
    for i in 1:length(pos)
        s.pos[i] .= pos[i]
    end

    for i in 2:s.set.segments+1
        state_y0[i-1] .= pos[i]  # Initial state vector
        yd0[i-1]      .= vel[i]  # Initial state vector derivative
    end

    for i in 2:s.set.segments+1
        state_y0[s.set.segments+i-1] .= vel[i]  # Initial state vector
        yd0[s.set.segments+i-1]      .= acc[i]  # Initial state vector derivative
    end
    set_v_wind_ground(s, pos[s.set.segments+1][3])
    s.l_tether = s.set.l_tether
    set_v_reel_out(s, s.set.v_reel_out, 0.0)

    state_y0, yd0
end

# same as above, but returns a tuple of two one dimensional arrays
function init(s::KPS3, X=zeros(2 * (s.set.segments)); old=false, delta = 0.0)
    res1_, res2_ = init_inner(s, X; old=old, delta = delta)
    res1, res2  = reduce(vcat, res1_), reduce(vcat, res2_)
    MVector{6*(s.set.segments), Float64}(res1), MVector{6*(s.set.segments), Float64}(res2)
end

"""
    spring_forces(s::AKM)

Return an array of the scalar spring forces of all tether segements.
"""
function spring_forces(s::KPS3)
    forces = zeros(SimFloat, s.set.segments)
    for i in 1:s.set.segments
        forces[i] =  s.c_spring * s.stiffness_factor * (norm(s.pos[i+1] - s.pos[i]) - s.segment_length)
    end
    forces
end

function find_steady_state_inner(s::KPS3, X, prn=false; delta=0.0)
    res = zeros(MVector{6*s.set.segments, SimFloat})

    # helper function for the steady state finder
    function test_initial_condition!(F, x::Vector)
        y0, yd0 = init(s, x, delta=delta)
        residual!(res, yd0, y0, s, 0.0)
        for i in 1:s.set.segments
            F[i]                = res[1 + 3*(i-1) + 3*s.set.segments]
            F[i+s.set.segments] = res[3 + 3*(i-1) + 3*s.set.segments]
        end
        return nothing 
    end

    if prn println("\nStarted function test_nlsolve...") end
    results = nlsolve(test_initial_condition!, X, xtol=1e-6, ftol=1e-6, autoscale=true, iterations=1000)
    if prn println("\nresult: $results") end
    results.zero
 end

"""
    find_steady_state(s::KPS3, prn=false, delta = 0.0, stiffness_factor=0.035)

Find an initial equilibrium, based on the inital parameters
`l_tether`, elevation and `v_reel_out`.
"""
function find_steady_state(s::KPS3; prn=false, delta = 0.0, stiffness_factor=0.035)
    zero = zeros(SimFloat, 2*s.set.segments)
    s.stiffness_factor=stiffness_factor
    zero = find_steady_state_inner(s, zero, prn, delta=delta)
    s.stiffness_factor=1.0
    zero = find_steady_state_inner(s, zero, prn, delta=delta)
    init(s, zero; delta=delta)
end