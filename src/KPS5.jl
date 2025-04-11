#= MIT License

Copyright (c) 2020, 2021, 2022, 2024 Uwe Fechner

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

This model implements a 3D mass-spring system with reel-out. It uses six tether segments (the number can be
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 3 aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on three of the four kite point masses. 

Four point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#

# Array of connections of bridlepoints.
# First point, second point, unstressed length.
const SPRINGS_INPUT_5P = [0.    1.  150.
                          1.    2.   -1. # s1, p7, p8
                          4.    2.   -1. # s2, p10, p8                        
                          4.    5.   -1. # s3, p10, p11
                          3.    4.   -1. # s4, p9, p10
                          5.    1.   -1. # s5, p11, p7
                          4.    1.   -1. # s6, p10, p7
                          3.    5.   -1. # s7, p9, p11
                          5.    2.   -1. # s8, p11, p8
                          2.    3.   -1.] # s9, p8, p9

# KCU = p7, A = p8, B = p9, C = p10, D = p11

"""
    mutable struct KPS5{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a 4 point kite model. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- T: Vector type, e.g. MVector{3, SimFloat}
- P: number of points of the system, segments+1
- Q: number of springs in the system, P-1
- SP: struct type, describing a spring
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPS5{S, T, P, Q, SP} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings
    "Reference to the KCU model (Kite Control Unit as implemented in the package KitePodModels"
    kcu::KCU
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "Reference to winch model as implemented in the package WinchModels"
    wm::AbstractWinchModel
    "Iterations, number of calls to the function residual!"
    iter:: Int64 = 0
    "Function for calculation the lift coefficent, using a spline based on the provided value pairs."
    calc_cl::Spline1D
    "Function for calculation the drag coefficent, using a spline based on the provided value pairs."
    calc_cd::Spline1D
    "wind vector at the height of the kite" 
    v_wind::T =           zeros(S, 3)
    "wind vector at reference height" 
    v_wind_gnd::T =       zeros(S, 3)
    "wind vector used for the calculation of the tether drag"
    v_wind_tether::T =    zeros(S, 3)
    "apparent wind vector at the kite"
    v_apparent::T =       zeros(S, 3)
    "bridle_factor = set.l_bridle/bridle_length(set)"
    bridle_factor::S = 1.0
    "side lift coefficient, the difference of the left and right lift coefficients"
    side_cl::S = 0.0
    "drag force of kite and bridle; output of calc_aero_forces!"
    drag_force::T =       zeros(S, 3)
    "side_force acting on the kite"
    side_force::T =       zeros(S, 3)
    "max_steering angle in radian"
    ks::S =               0.0
    "lift force of the kite; output of calc_aero_forces!"
    lift_force::T =       zeros(S, 3)    
    "spring force of the current tether segment, output of calc_particle_forces!"
    spring_force::T =     zeros(S, 3)
    "last winch force"
    last_force::T =       zeros(S, 3)
    "a copy of the residual one (pos,vel) for debugging and unit tests"    
    res1::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the residual two (vel,acc) for debugging and unit tests"
    res2::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the actual positions as output for the user"
    pos::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "a copy of the actual velocities as output for the user"
    vel::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "velocity vector of the kite"
    vel_kite::T =          zeros(S, 3)
    "unstressed segment length [m]"
    segment_length::S =           0.0
    "lift coefficient of the kite, depending on the angle of attack"
    param_cl::S =         0.2
    "drag coefficient of the kite, depending on the angle of attack"
    param_cd::S =         1.0
    "azimuth angle in radian; initial value is zero"
    psi::S =              zero(S)
    "depower angle [deg]"
    alpha_depower::S =     0.0
    "pitch angle [rad]"
    pitch::S =            0.0
    "pitch rate [rad/s]"
    pitch_rate::S =       0.0
    "aoa at particle B"
    alpha_2::S =           0.0
    "aoa at particle B, corrected formula"
    alpha_2b::S =          0.0
    "aoa at particle C"
    alpha_3::S =           0.0
    alpha_3b::S =          0.0
    "aoa at particle D"
    alpha_4::S =           0.0
    alpha_4b::S =          0.0
    "relative start time of the current time interval"
    t_0::S =               0.0
    "reel out speed of the winch"
    v_reel_out::S =        0.0
    "reel out speed at the last time step"
    last_v_reel_out::S =   0.0
    "unstretched tether length"
    l_tether::S =          0.0
    "air density at the height of the kite"
    rho::S =               0.0
    "actual relative depower setting,  must be between    0 .. 1.0"
    depower::S =           0.0
    "actual relative steering setting, must be between -1.0 .. 1.0"
    steering::S =          0.0
    "steering after the kcu, before applying offset and depower sensitivity, -1.0 .. 1.0"
    kcu_steering::S =      0.0
    "multiplier for the stiffniss of tether and bridle"
    stiffness_factor::S =  1.0
    "initial masses of the point masses"
    initial_masses::MVector{P, S} = ones(P)
    "current masses, depending on the total tether length"
    masses::MVector{P, S}         = zeros(P)
    "vector of the springs, defined as struct"
    springs::MVector{Q, SP}       = zeros(SP, Q)
    "vector of the forces, acting on the particles"
    forces::SVector{P, KVec3} = zeros(SVector{P, KVec3})
    "synchronous speed of the motor/ generator"
    sync_speed::Union{S, Nothing} =        0.0
    "set_torque of the motor/generator"
    set_torque::Union{S, Nothing} = nothing
    "set value of the force at the winch, for logging only"
    set_force::Union{S, Nothing} = nothing
    "set value of the bearing angle in radians, for logging only"
    bearing::Union{S, Nothing} = nothing
    "coordinates of the attractor point [azimuth, elevation] in radian, for logging only"
    attractor::Union{SVector{2, S}, Nothing} = nothing
    "x vector of kite reference frame"
    x::T =                 zeros(S, 3)
    "y vector of kite reference frame"
    y::T =                 zeros(S, 3)
    "z vector of kite reference frame"
    z::T =                 zeros(S, 3)
end

"""
    clear!(s::KPS5)

Initialize the kite power model.
"""
function clear!(s::KPS5)
    s.t_0 = 0.0                              # relative start time of the current time interval
    s.v_reel_out = s.set.v_reel_out
    s.last_v_reel_out = s.set.v_reel_out
    s.sync_speed = s.set.v_reel_out
    s.v_wind_gnd    .= [s.set.v_wind, 0, 0]    # wind vector at reference height
    s.v_wind_tether .= [s.set.v_wind, 0, 0]
    s.v_apparent    .= [s.set.v_wind, 0, 0]
    height = sin(deg2rad(s.set.elevation)) * (s.set.l_tether)
    s.v_wind .= s.v_wind_gnd * calc_wind_factor(s.am, height)

    s.l_tether = s.set.l_tether
    s.segment_length = s.l_tether / s.set.segments
    init_masses!(s)
    init_springs!(s)
    for i in 1:s.set.segments + KiteModels.KITE_PARTICLES + 1 
        s.forces[i] .= zeros(3)
    end
    s.drag_force .= [0.0, 0, 0]
    s.lift_force .= [0.0, 0, 0]
    s.side_force .= [0.0, 0, 0]
    s.rho = s.set.rho_0
    s.bridle_factor = s.set.l_bridle / bridle_length(s.set)
    s.ks = deg2rad(s.set.max_steering) 
    s.kcu.depower = s.set.depower/100.0
    s.kcu.set_depower = s.kcu.depower
    roll, pitch, yaw = orient_euler(s)
    s.pitch = pitch
    s.pitch_rate = 0.0
    KiteModels.set_depower_steering!(s, get_depower(s.kcu), get_steering(s.kcu))
end

function KPS5(kcu::KCU)
    if kcu.set.winch_model == "AsyncMachine"
        wm = AsyncMachine(kcu.set)
    elseif kcu.set.winch_model == "TorqueControlledMachine"
        wm = TorqueControlledMachine(kcu.set)
    end
    # wm.last_set_speed = kcu.set.v_reel_out
    s = KPS5{SimFloat, KVec3, kcu.set.segments+KITE_PARTICLES+1, kcu.set.segments+KITE_SPRINGS, SP}(set=kcu.set, 
             kcu=kcu, wm=wm, calc_cl = Spline1D(kcu.set.alpha_cl, kcu.set.cl_list), 
             calc_cd=Spline1D(kcu.set.alpha_cd, kcu.set.cd_list) )    
    clear!(s)
    return s
end

""" 
    calc_particle_forces!(s::KPS5, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)

Calculate the drag force of the tether segment, defined by the parameters pos1, pos2, vel1 and vel2
and distribute it equally on the two particles, that are attached to the segment.
The result is stored in the array s.forces. 
"""
@inline function calc_particle_forces!(s::KPS5, pos1, pos2, vel1, vel2, spring, segments, d_tether, rho, i)
    l_0 = spring.length # Unstressed length
    k = spring.c_spring * s.stiffness_factor  # Spring constant
    c = spring.damping                        # Damping coefficient    
    segment = pos1 - pos2
    rel_vel = vel1 - vel2
    av_vel = 0.5 * (vel1 + vel2)
    norm1 = norm(segment)
    unit_vector = segment / norm1

    k1 = s.set.rel_compr_stiffness * k # compression stiffness kite springs
    k2 = 0.1 * k                       # compression stiffness tether springs
    c1 = s.set.rel_damping * c         # damping kite springs
    spring_vel   = unit_vector ⋅ rel_vel
    if (norm1 - l_0) > 0.0
        if i > segments  # kite springs
             s.spring_force .= (k *  (norm1 - l_0) + (c1 * spring_vel)) * unit_vector 
        else
             s.spring_force .= (k *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
        end
    elseif i > segments # kite spring
        s.spring_force .= (k1 *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
    else
        s.spring_force .= (k2 *  (norm1 - l_0) + (c * spring_vel)) * unit_vector
    end

    s.v_apparent .= s.v_wind_tether - av_vel
    v_app_kcu = s.v_wind_tether - vel2
    if s.set.version == 1
        area = norm1 * d_tether
    else
        if i > segments
            area = norm1 * s.set.d_line * 0.001 * s.bridle_factor # 6.0 = A_real/A_simulated
        else
            area = norm1 * d_tether
        end
    end

    v_app_perp = s.v_apparent - s.v_apparent ⋅ unit_vector * unit_vector
    half_drag_force = (-0.25 * rho * s.set.cd_tether * norm(v_app_perp) * area) * v_app_perp 
    if i == segments
        v_app_perp_kcu = v_app_kcu - v_app_kcu ⋅ unit_vector * unit_vector
        kcu_area = π * (s.set.kcu_diameter/2)^2
        kcu_drag_force = (-0.25 * rho * s.set.cd_kcu * norm(v_app_perp_kcu) * kcu_area) * v_app_perp_kcu
        @inbounds s.forces[spring.p2] .+= kcu_drag_force
    end

    @inbounds s.forces[spring.p1] .+= half_drag_force + s.spring_force
    @inbounds s.forces[spring.p2] .+= half_drag_force - s.spring_force
    if i == 1 s.last_force .= s.forces[spring.p1] end
    nothing
end

"""
    calc_aero_forces!(s::KPS5, pos, vel, rho, alpha_depower, rel_steering)

Calculates the aerodynamic forces acting on the kite particles.

Parameters:
- pos:              vector of the particle positions
- vel:              vector of the particle velocities
- rho:              air density [kg/m^3]
- rel_depower:      value between  0.0 and  1.0
- alpha_depower:    depower angle [degrees]
- rel_steering:     value between -1.0 and +1.0

Updates the vector s.forces of the first parameter.
"""
function calc_aero_forces!(s::KPS5, pos, vel, rho, alpha_depower, rel_steering)
end


# =================== getter functions ====================================================

"""
    calc_height(s::KPS5)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPS5)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPS5)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPS5)
    s.pos[end-2]
end

"""
    kite_ref_frame(s::KPS5; one_point=false)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
"""
function kite_ref_frame(s::KPS5; one_point=false)
    if one_point
        c = s.z
        y = normalize(s.v_apparent × c)
        x = normalize(y × c)
        return x, y, c
    else
        return s.x, s.y, s.z
    end
end

"""
    winch_force(s::KPS5)

Return the absolute value of the force at the winch as calculated during the last timestep. 
"""
function winch_force(s::KPS5) norm(s.last_force) end

"""
    cl_cd(s::KPS5)

Calculate the lift and drag coefficients of the kite, based on the current angles of attack.
"""
function cl_cd(s::KPS5)
    rel_side_area = s.set.rel_side_area/100.0  # defined in percent
    K = 1 - rel_side_area                      # correction factor for the drag
    if s.set.version == 3
        drag_corr = 1.0
    else
        drag_corr = DRAG_CORR
    end
    CL2, CD2 = s.calc_cl(s.alpha_2), drag_corr * s.calc_cd(s.alpha_2)
    CL3, CD3 = s.calc_cl(s.alpha_3), drag_corr * s.calc_cd(s.alpha_3)
    CL4, CD4 = s.calc_cl(s.alpha_4), drag_corr * s.calc_cd(s.alpha_4)
    if s.set.version == 3
        return CL2, CD2
    else
        return CL2, K*(CD2+CD3+CD4)
    end
end

# ==================== end of getter functions ================================================

function spring_forces(s::KPS5; prn=true)
    forces = zeros(SimFloat, s.set.segments+KITE_SPRINGS)
    for i in 1:s.set.segments
        forces[i] =  s.springs[i].c_spring * (norm(s.pos[i+1] - s.pos[i]) - s.segment_length) * s.stiffness_factor
        if forces[i] > s.set.max_force && prn
            println("Tether raptures for segment $i !")
        end
    end
    for i in 1:KITE_SPRINGS
        p1 = s.springs[i+s.set.segments].p1  # First point nr.
        p2 = s.springs[i+s.set.segments].p2  # Second point nr.
        pos1, pos2 = s.pos[p1], s.pos[p2]
        spring = s.springs[i+s.set.segments]
        l_0 = spring.length # Unstressed length
        k = spring.c_spring * s.stiffness_factor       # Spring constant 
        segment = pos1 - pos2
        norm1 = norm(segment)
        k1 = 0.25 * k # compression stiffness kite segments
        if (norm1 - l_0) > 0.0
            spring_force = k *  (norm1 - l_0) 
        else 
            spring_force = k1 *  (norm1 - l_0)
        end
        forces[i+s.set.segments] = spring_force
        if norm(s.spring_force) > 4000.0
            println("Bridle brakes for spring $i !")
        end
    end
    forces
end

"""
    find_steady_state!(s::KPS5; prn=false, delta = 0.01, stiffness_factor=0.035, upwind_dir=-pi/2))

Find an initial equilibrium, based on the initial parameters
`l_tether`, elevation and `v_reel_out`.
"""
function find_steady_state!(s::KPS5; prn=false, delta = 0.01, stiffness_factor=0.035, upwind_dir=-pi/2)
    set_v_wind_ground!(s, calc_height(s), s.set.v_wind; upwind_dir=-pi/2)
    s.stiffness_factor = stiffness_factor
    res = zeros(MVector{6*(s.set.segments+KITE_PARTICLES)+2, SimFloat})
    iter = 0

    # helper function for the steady state finder
    function test_initial_condition!(F, x::Vector)
        x1 = copy(x)
        y0, yd0 = init(s, x1; delta)
        try
            residual!(res, yd0, y0, s, 0.0)
        catch e
            println("Warning in test_initial_condition!")
        end
        for i in 1:s.set.segments+KITE_PARTICLES-1
            if i != s.set.segments+KITE_PARTICLES-1
                j = i
            else
                j = i + 1
            end
            # copy the x-component of the residual res2 (acceleration)
            F[i]                               = res[1 + 3*(j-1) + 3*(s.set.segments+KITE_PARTICLES)]
            # copy the z-component of the residual res2
            F[i+s.set.segments+KITE_PARTICLES] = res[3 + 3*(j-1) + 3*(s.set.segments+KITE_PARTICLES)]
        end
        # copy the acceleration of point KCU in x direction
        i = s.set.segments+1
        F[end-1]                               = res[1 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)] 
        # copy the acceleration of point C in y direction
        i = s.set.segments+3 
        x = res[1 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)]
        y = res[2 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)]
        F[end]                                 = res[2 + 3*(i-1) + 3*(s.set.segments+KITE_PARTICLES)] 
        iter += 1
        return nothing 
    end
    if prn println("\nStarted function test_nlsolve...") end
    X00 = zeros(SimFloat, 2*(s.set.segments+KITE_PARTICLES-1)+2)
    results = nlsolve(test_initial_condition!, X00, autoscale=true, xtol=4e-7, ftol=4e-7, iterations=s.set.max_iter)
    if prn println("\nresult: $results") end
    y0, yd0 = init(s, results.zero; upwind_dir)
    set_v_wind_ground!(s, calc_height(s), s.set.v_wind; upwind_dir)
    residual!(res, yd0, y0, s, 0.0)
    y0, yd0
end
