#= MIT License

Copyright (c) 2024 Uwe Fechner and Bart van de Lint

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
configured in the file data/settings.yaml). The kite is modelled using 4 point masses and 2n aerodynamic 
surfaces. The spring constant and the damping decrease with the segment length. The aerodynamic kite forces
are acting on the inertial kite point. 

Four point kite model, included from KiteModels.jl.

Scientific background: http://arxiv.org/abs/1406.6218 =#

# Array of connections of bridlepoints.
# First point, second point, unstressed length.
const SPRINGS_INPUT_3L = [1.      4.  -1. # s1: E, A
                        1.      2.  -1. # s2, E, C
                        1.      3.  -1. # s3, E, D
                        2.      3.  -1. # s4, C, D
                        2.      4.  -1. # s5, C, A
                        3.      4.  -1. # s6, D, A
                        ]
# E = 1, C = 2, D = 3, A = 4
# E = segments*3+1, C = segments*3+2, D = segments*3+3, A = segments*3+4

# TODO: add line multiplier: multiple lines doing same thing
const KITE_SPRINGS_3L = 6
const KITE_PARTICLES_3L = 4
const MeasureFloat = Float32

@with_kw mutable struct Measurement
    set_values::MVector{3, MeasureFloat}    = [-1., -1., -50.]
    tether_length::MVector{3, MeasureFloat} = [51., 51., 49.]
    tether_vel::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_acc::MVector{3, MeasureFloat}    = zeros(MeasureFloat, 3)
    tether_force::MVector{3, MeasureFloat}  = [3., 3., 540.]
    "elevation and azimuth in spherical coordinate system with columns (left, right) and rows (elevation, azimuth)"
    sphere_pos::Matrix{MeasureFloat}            = deg2rad.([89.0 89.0; 1.0 -1.0])
    sphere_vel::Matrix{MeasureFloat}            = zeros(MeasureFloat, 2, 2)
    sphere_acc::Matrix{MeasureFloat}            = zeros(MeasureFloat, 2, 2)
    "positive azimuth wind direction in right-handed ENU frame relative to east / x-axis"
    wind_dir_gnd::MeasureFloat                  = zero(MeasureFloat)
end

function Base.getproperty(m::Measurement, val::Symbol)
    if val === :elevation
        sphere_pos = getfield(m, :sphere_pos)
        return 0.5(sphere_pos[1, 1] + sphere_pos[1, 2])
    elseif val === :azimuth
        sphere_pos = getfield(m, :sphere_pos)
        return 0.5(sphere_pos[2, 1] + sphere_pos[2, 2])
    else
        return getfield(m, val)
    end
end

@enum SegmentType begin
    POWER
    STEERING
    BRIDLE
end

@enum DynamicsType begin
    DYNAMIC
    STATIC
    KITE
    WINCH
end

"""
A normal freely moving tether point
"""
mutable struct Point
    idx::Int16
    pos_b::KVec3 # pos relative to kite COM in body frame
    pos_w::KVec3 # pos in world frame
    type::DynamicsType
end
function Point(idx, pos_b, type)
    Point(idx, pos_b, copy(pos_b), type)
end

"""
Set of bridle lines that share the same twist angle and trailing edge angle
"""
struct KitePointGroup
    idx::Int16
    points::Vector{Int16}
    y_lim::Tuple{SimFloat, SimFloat}
    fixed_index::Int16 # point which the group rotates around under kite deformation
    chord::KVec3 # chord vector in body frame which the group rotates around under kite deformation
    y_airf::KVec3 # spanwise vector in local panel frame which the group rotates around under kite deformation
    type::DynamicsType
end

"""
A segment from one point index to another point index
"""
mutable struct Segment
    idx::Int16
    points::Tuple{Int16, Int16}
    type::SegmentType
    l0::SimFloat
    diameter::SimFloat
end
function Segment(idx, points, type)
    Segment(idx, points, type, zero(SimFloat), zero(SimFloat))
end
function Segment(idx, points, type, l0)
    Segment(idx, points, type, l0, zero(SimFloat))
end

"""
A pulley described by two segments with the common point of the segments being the pulley
"""
mutable struct Pulley
    idx::Int16
    segments::Tuple{Int16, Int16}
    type::DynamicsType
    sum_length::SimFloat
    function Pulley(idx, segments, type)
        new(idx, segments, type, zero(SimFloat))
    end
end

"""
A set of segments making a flexible tether. The winch point should only be part of one segment.
"""
struct Tether
    idx::Int16
    segments::Vector{Int16}
    winch_point::Int16
end

"""
A set of tethers or just one tether connected to a winch
"""
mutable struct Winch
    idx::Int16
    model::AbstractWinchModel
    tethers::Vector{Int16}
    tether_length::Float64
    function Winch(idx, model, tethers)
        new(idx, model, tethers, zero(Float64))
    end
end

struct PointMassSystem
    points::Vector{Point}
    groups::Vector{KitePointGroup}
    segments::Vector{Segment}
    pulleys::Vector{Pulley}
    tethers::Vector{Tether}
    winches::Vector{Winch}
    function PointMassSystem(points, groups, segments, pulleys, tethers, winches)
        for (i, point) in enumerate(points)
            @assert point.idx == i
        end
        for (i, group) in enumerate(groups)
            @assert group.idx == i
        end
        for (i, segment) in enumerate(segments)
            @assert segment.idx == i
        end
        for (i, pulley) in enumerate(pulleys)
            @assert pulley.idx == i
        end
        for (i, tether) in enumerate(tethers)
            @assert tether.idx == i
        end
        for (i, winch) in enumerate(winches)
            @assert winch.idx == i
        end
        new(points, groups, segments, pulleys, tethers, winches)
    end
end

"""
    mutable struct KPSQ{S, T, P, Q, SP} <: AbstractKiteModel

State of the kite power system, using a quaternion kite model and three steering lines to the ground. Parameters:
- S: Scalar type, e.g. SimFloat
  In the documentation mentioned as Any, but when used in this module it is always SimFloat and not Any.
- V: Vector type, e.g. KVec3
- P: number of tether points of the system, 3segments+3
Normally a user of this package will not have to access any of the members of this type directly,
use the input and output functions instead.

$(TYPEDFIELDS)
"""
@with_kw mutable struct KPSQ{S, V, P} <: AbstractKiteModel
    "Reference to the settings struct"
    set::Settings
    "Reference to the geometric wing model"
    wing::VortexStepMethod.RamAirWing
    "Reference to the aerodynamic wing model"
    aero::VortexStepMethod.BodyAerodynamics
    "Reference to the VSM aerodynamics solver"
    vsm_solver::VortexStepMethod.Solver
    "Reference to the point mass system with points, segments, pulleys and tethers"
    point_system::PointMassSystem = PointMassSystem(Point[], KitePointGroup[], Segment[], Pulley[], Tether[], Winch[])
    "The last initial elevation"
    last_init_elevation::S     = 0.0
    "The last initial tether length"
    last_init_tether_length::S = 0.0
    "Reference to the last settings hash"
    last_set_hash::UInt64   = 0
    "Reference to the last measurement hash"
    last_measure_hash::UInt64 = 0
    "Reference to the atmospheric model as implemented in the package AtmosphericModels"
    am::AtmosphericModel = AtmosphericModel()
    "tether positions"
    pos::Matrix{S} = zeros(S, 3, P)
    "unstressed segment lengths of the three tethers [m]"
    segment_lengths::V =           zeros(S, 3)
    "relative start time of the current time interval"
    t_0::S =               0.0
    "unstretched tether length"
    tether_lengths::V =          zeros(S, 3)
    "air density at the height of the kite"
    rho::S =               0.0
    "tether masses"
    masses::V         = zeros(S, P)
    "unit spring coefficient"
    c_spring::V = zeros(S, 3)
    "unit damping coefficient"
    damping::V = zeros(S, 3)
    "whether or not to use torque control instead of speed control"
    torque_control::Bool = false
    "x vector of kite reference frame"
    e_x::V =                 zeros(S, 3)
    "y vector of kite reference frame"
    e_y::V =                 zeros(S, 3)
    "z vector of kite reference frame"
    e_z::V =                 zeros(S, 3)
    "Simplified system of the mtk model"
    simple_sys::Union{ModelingToolkit.ODESystem, Nothing} = nothing
    "Velocity of the kite"
    vel_kite::V =           zeros(S, 3)
    "Initial torque or speed set values"
    init_set_values::V =    zeros(S, 3)
    "Smooth sign constant"
    ϵ::S =      0.0
    "Relative damping for flaps"
    flap_damping::S     = 0.75
    "Measured data points used to create an initial state"
    measure::Measurement = Measurement()
    "Makes autodiff faster"
    prep::Union{Nothing, Any} = nothing
    "Buffer for jacobian y values (vectorized velocities)"
    y_buffer::V = zeros(S, P)
    "Buffer for jacobian x values (angle, distance)"
    x_buffer::V = zeros(S, 2)
    "Inertia around kite x y and z axis of the body frame"
    I_b::V = zeros(S, 3)
    "Damping of the kite rotation"
    orient_damping::S = zero(S)
    "Initialization values for kite state"
    u0map::Union{Vector{Pair{Num, S}}, Nothing} = nothing
    "Initialization values for kite parameters"
    p0map::Union{Vector{Pair{Num, S}}, Nothing} = nothing
    "Distance of the kite com from winch"
    distance::S = zero(S)
    "X coordinate on normalized 2d foil of bridle attachments"
    bridle_fracs::V = [0.088, 0.31, 0.58, 0.93]
    crease_frac::S = 0.82
    "The top bridle points that are not on the kite, in CAD frame"
    top_bridle_points::Vector{V} = [[0.290199, 0.784697, -2.61305], [0.392683, 0.785271, -2.61201], [0.498202, 0.786175, -2.62148], [0.535543, 0.786175, -2.62148]]
    "Tether diameter of tethers in bridle system [mm]"
    bridle_tether_diameter::SimFloat = 2.
    "Tether diameter of the power tethers [mm]"
    power_tether_diameter::SimFloat = 2.
    "Tether diameter of the steering tethers [mm]"
    steering_tether_diameter::SimFloat = 1.

    set_set_values::Function       = () -> nothing
    set_measure::Function          = () -> nothing
    set_coefficients::Function     = () -> nothing

    get_state::Function            = () -> nothing
    get_twist::Function            = () -> nothing
    get_va_body::Function          = () -> nothing

    prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    init_prob::Union{OrdinaryDiffEqCore.ODEProblem, Nothing} = nothing
    integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Sundials.CVODEIntegrator, Nothing} = nothing
    init_integrator::Union{OrdinaryDiffEqCore.ODEIntegrator, Sundials.CVODEIntegrator, Nothing} = nothing
end

"""
    clear!(s::KPSQ)

Initialize the kite power model.
"""
function clear!(s::KPSQ)
    # P = 3s.set.segments + 3
    # S = eltype(s.pos)
    # s.pos = zeros(S, 3, P)
    # s.masses = zeros(S, P)
    # s.expected_tether_pos_vel_buffer = zeros(S, 2, 3, (P ÷ 3))
    # s.J_buffer = zeros(S, P, 2)
    # s.y_buffer = zeros(S, P)
    # s.t_0 = 0.0                              # relative start time of the current time interval
    # s.e_x .= 0.0
    # s.e_y .= 0.0
    # s.e_z .= 0.0
    # s.tether_lengths .= [s.set.l_tether for _ in 1:3]
    # s.γ_l = π/2 - s.set.min_steering_line_distance/(2*s.set.radius)
    # s.segment_lengths .= s.tether_lengths ./ s.set.segments
    # s.i_A = s.set.segments*3+1
    # s.i_B = s.set.segments*3+2
    # s.i_C = s.set.segments*3+3
    # s.rho = s.set.rho_0
    # c_spring = s.set.e_tether * (s.set.d_tether/2000.0)^2 * pi
    # s.c_spring .= [c_spring, c_spring, 2c_spring]
    # s.damping .= (s.set.damping / s.set.c_spring) * s.c_spring
    # s.kite_length = function (γ)
    #     if γ < π/2
    #         return s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (γ - s.γ_l) / (π/2 - s.γ_l)
    #     else
    #         return s.set.tip_length + (s.set.middle_length-s.set.tip_length) * (π - s.γ_l - γ) / (π/2 - s.γ_l)
    #     end
    # end
    # init_masses!(s)

    # width, radius, tip_length, middle_length = s.set.width, s.set.radius, s.set.tip_length, s.set.middle_length
    # s.γ_l = pi/2 - width/2/radius
    # s.γ_D = s.γ_l + width*(-2*tip_length + sqrt(2*middle_length^2 + 2*tip_length^2)) /
    #     (4*(middle_length - tip_length)) / radius
    # s.kite_length_D = tip_length + (middle_length-tip_length) * (s.γ_D - s.γ_l) / (π/2 - s.γ_l)

    # calc_inertia!(s, s.wing)
    # calc_pos_principal!(s)
    nothing
end

function KPSQ(set::Settings, wing::RamAirWing, aero::BodyAerodynamics, vsm_solver::VortexStepMethod.Solver)
    if set.winch_model == "TorqueControlledMachine"
        s = KPSQ{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            ; set, wing, aero, vsm_solver
            )
        s.torque_control = true
    else
        s = KPSQ{SimFloat, Vector{SimFloat}, 3*(set.segments + 1)}(
            ; set, wing, aero, vsm_solver
            )
        s.torque_control = false
    end
    clear!(s)
    return s
end

function calc_kite_ref_frame!(s::KPSQ, E, C, D)
    P_c = 0.5 .* (C+D)
    s.e_y .= normalize(C - D)
    s.e_z .= normalize(E - P_c)
    s.e_x .= cross(s.e_y, s.e_z)
    return nothing
end

function calc_tether_elevation(s::KPSQ)
    KiteUtils.calc_elevation(s.pos[6])
end

function calc_tether_azimuth(s::KPSQ)
    KiteUtils.azimuth_east(s.pos[6])
end

function update_sys_state!(ss::SysState, s::KPSQ, zoom=1.0)
    ss.time = s.t_0
    pos, acc, Q_p_w, elevation, azimuth, e_x, tether_vel, twist, kite_vel = s.get_state(s.integrator)
    P = length(s.point_system.points)
    for i in 1:P
        ss.X[i] = pos[1, i] * zoom
        ss.Y[i] = pos[2, i] * zoom
        ss.Z[i] = pos[3, i] * zoom
    end
    # TODO
    # ss.kite_acc      .= kite_acc(s)
    # ss.left_tether_vel = tether_vel[1]
    # ss.right_tether_vel = tether_vel[2]
    ss.acc = norm(acc)
    ss.orient .= Q_p_w
    ss.elevation = elevation
    ss.azimuth = azimuth
    ss.force = zero(SimFloat)
    ss.heading = calc_heading_y(e_x)
    ss.course = calc_course(1.0, 1.0) # TODO: implement
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = mean(tether_vel)
    ss.depower = rad2deg(mean(twist))
    ss.steering = rad2deg(twist[end] - twist[1])
    ss.vel_kite .= kite_vel
    nothing
end

function SysState(s::KPSQ, zoom=1.0) # TODO: add left and right lines, stop using getters and setters
    isnothing(s.integrator) && throw(ArgumentError("run init!(s) first"))
    pos, acc, Q_p_w, elevation, azimuth, e_x, tether_vel, twist, kite_vel = s.get_state(s.integrator)
    P = length(s.point_system.points)
    X = zeros(MVector{P, MyFloat})
    Y = zeros(MVector{P, MyFloat})
    Z = zeros(MVector{P, MyFloat})
    for i in 1:P
        X[i] = pos[1, i] * zoom
        Y[i] = pos[2, i] * zoom
        Z[i] = pos[3, i] * zoom
    end
    
    orient = MVector{4, Float32}(Q_p_w) # TODO: add Q_b_w
    # forces = s.get_tether_force() # TODO: add tether force
    forces = zeros(3)
    heading = calc_heading_y(e_x)
    course = calc_course(1.0, 1.0) # TODO: implement
    t_sim = 0
    depower = rad2deg(mean(twist))
    steering = rad2deg(twist[end] - twist[1])
    ss = SysState{P}()
    ss.time = s.t_0
    ss.t_sim = t_sim
    ss.orient .= orient
    ss.elevation = elevation
    ss.azimuth = azimuth
    ss.l_tether = s.tether_lengths[3]
    ss.v_reelout = tether_vel[3]
    ss.force = forces[3]
    ss.depower = depower
    ss.steering = steering
    ss.heading = heading
    ss.course = course
    ss.vel_kite .= kite_vel
    ss.X = X
    ss.Y = Y
    ss.Z = Z
    ss
end


# function calc_heading(e_x, pos_kite)
#     # turn s.e_x by -azimuth around global z-axis and then by elevation around global y-axis
#     vec = rotate_around_y(rotate_around_z(e_x, -KiteUtils.azimuth_east(pos_kite)), -KiteUtils.calc_elevation(pos_kite))
#     heading = atan(-vec[2], vec[3])
#     return heading
# end
# @register_symbolic calc_heading(e_x, pos_kite)

function calc_heading_y(e_x)
    return atan(-e_x[2]/-e_x[1])
end
@register_symbolic calc_heading_y(e_x)

function calc_elevation(s::KPSQ) return s.get_elevation() end
function calc_azimuth(s::KPSQ) return s.get_azimuth() end

"""
course where straight up is zero, clockwise is positive
"""
function calc_course(elevation_vel, azimuth_vel)
    course = atan(-azimuth_vel, elevation_vel)
    return course
end

"""
Initialises the integrator of the model.

Parameters:
- s:     an instance of a 3 line kite model
- prn: if set to true, print the detailed solver results
- torque_control: wether or not to use torque control

Returns:
Nothing.
"""
function init_sim!(s::KPSQ; prn=false, torque_control=s.torque_control, 
        init_set_values=s.init_set_values, ϵ=s.ϵ, flap_damping=s.flap_damping, 
        force_new_sys=false, force_new_pos=false, init=false)
    dt = SimFloat(1/s.set.sample_freq)
    tspan   = (0.0, dt) 
    solver = Rodas5P(autodiff=AutoFiniteDiff())
    set_hash = struct_hash(s.set)
    measure_hash = struct_hash(s.measure)
    new_pos = s.last_measure_hash != measure_hash || force_new_pos
    new_sys = force_new_sys ||
                    isnothing(s.prob) || 
                    s.torque_control != torque_control || 
                    s.last_set_hash != set_hash ||
                    s.ϵ != ϵ ||
                    s.flap_damping != flap_damping
    s.torque_control = torque_control
    s.ϵ = ϵ
    s.flap_damping = flap_damping

    if new_sys
        if prn; println("initializing with new model and new pos"); end
        clear!(s)
        sys, defaults, guesses = model!(s; init)
        # for u in unknowns(sys)
        #     println(u)
        # end
        @info "Creating the problem"
        @time s.prob = ODEProblem(sys, defaults, tspan; guesses)
        s.simple_sys = s.prob.f.sys
        s.integrator = OrdinaryDiffEqCore.init(s.prob, solver; dt, abstol=s.set.abs_tol, reltol=s.set.rel_tol, save_on=false)
        generate_getters!(s; init)
        if (!init) reinit!(s; new_sys) end
    elseif new_pos
        if prn; println("initializing with last model and new pos"); end
        if (!init) reinit!(s; new_sys) end
    else
        if prn; println("initializing with last model and last pos"); end
        OrdinaryDiffEqCore.reinit!(s.integrator, s.prob.u0)
    end

    s.last_init_elevation = s.set.elevation
    s.last_init_tether_length = s.set.l_tether
    s.last_set_hash = set_hash
    s.last_measure_hash = measure_hash
    s.init_set_values .= init_set_values
    return nothing
end

function generate_getters!(s; init=false)
    sys = s.prob.f.sys
    c = collect
    set_set_values = setp(sys, sys.set_values)
    set_measure = setp(sys, sys.measured_wind_dir_gnd)
    set_coefficients = setp(sys, [
        sys.moment_dist,
        sys.aero_kite_force_b,
        sys.aero_kite_moment_b
    ])
    get_state = getu(sys, 
        [c(sys.pos), c(sys.acc), c(sys.Q_p_w), sys.elevation, sys.azimuth, 
        c(sys.e_x), c(sys.tether_vel), c(sys.twist_angle), c(sys.kite_vel)]
    )
    get_twist = getu(sys, sys.twist_angle)
    get_va_body = getu(sys, sys.va_kite_b)

    s.set_set_values = (integ, val) -> set_set_values(integ, val)
    s.set_measure = (integ, val) -> set_measure(integ, val)
    s.set_coefficients = (integ, val) -> set_coefficients(integ, val)
    s.get_state = (integ) -> get_state(integ)
    s.get_twist = (integ) -> get_twist(integ)
    s.get_va_body = (integ) -> get_va_body(integ)
    nothing
end

function refine_twist!(s::KPSQ, twist_dist)
    twist_angles = s.get_twist(s.integrator)
    groups = s.point_system.groups
    panels = s.aero.panels
    group_idx = 1
    for (i, panel) in enumerate(panels)
        @show 
        if 0.5(panel.LE_point_2[2] + panel.LE_point_1[2]) < groups[group_idx].y_lim[2]
            group_idx += 1
        end
        (group_idx > length(groups)) && throw(ArgumentError("Panels and groups are not sorted in the same direction"))
        twist_dist[i] = twist_angles[group_idx]
    end
    @assert (group_idx == length(groups))

    window_size = length(panels) ÷ length(twist_angles)
    if length(panels) > window_size
        smoothed = copy(twist_dist)
        for i in (window_size÷2 + 1):(length(panels) - window_size÷2)
            smoothed[i] = mean(twist_dist[(i - window_size÷2):(i + window_size÷2)])
        end
        twist_dist .= smoothed
    end
    return nothing
end

function next_step!(s::KPSQ; set_values=nothing, measure::Union{Measurement, Nothing}=nothing, dt=1/s.set.sample_freq)
    if (!isnothing(set_values)) 
        s.set_set_values(s.integrator, set_values)
    end
    if (!isnothing(measure))
        s.set_measure(s.integrator, s.measure.wind_dir_gnd)
    end

    twist_distribution = zeros(length(s.aero.panels))
    refine_twist!(s, twist_distribution)
    VortexStepMethod.deform!(s.wing, twist_distribution, zeros(length(s.aero.panels)))
    VortexStepMethod.init!(s.aero)

    va_body = s.get_va_body(s.integrator)
    VortexStepMethod.set_va!(s.aero, va_body)
    VortexStepMethod.solve!(s.vsm_solver, s.aero; moment_frac=s.bridle_fracs[s.point_system.groups[1].fixed_index])
    if !any(isnan.(va_body)) &&
        !any(isnan.(s.vsm_solver.sol.gamma_distribution))

        s.set_coefficients(s.integrator, [
            s.vsm_solver.sol.moment_distribution,
            s.vsm_solver.sol.aero_force,
            s.vsm_solver.sol.aero_moments
        ])
    else
        @warn "Not converged"
    end
    
    s.t_0 = s.integrator.t
    OrdinaryDiffEqCore.step!(s.integrator, dt, true)
    if !successful_retcode(s.integrator.sol)
        println("Return code for solution: ", s.integrator.sol.retcode)
    end
    @assert successful_retcode(s.integrator.sol)
    s.integrator.t
end

"""
    unstretched_length(s::KPSQ)

Getter for the unstretched tether reel-out length (at zero force).
"""
function unstretched_length(s::KPSQ) s.tether_lengths[3] end


# =================== getter functions ====================================================

"""
    calc_height(s::KPSQ)

Determine the height of the topmost kite particle above ground.
"""
function calc_height(s::KPSQ)
    pos_kite(s)[3]
end

"""
    pos_kite(s::KPSQ)

Return the position of the kite (top particle).
"""
function pos_kite(s::KPSQ)
    s.kite_pos
end

"""
    kite_ref_frame(s::KPSQ; one_point=false)

Returns a tuple of the x, y, and z vectors of the kite reference frame.
The parameter one_point is not used in this model.
"""
function kite_ref_frame(s::KPSQ; one_point=false)
    s.get_e_x(), s.get_e_y(), s.get_e_z()
end

# ====================== helper functions ====================================

function struct_hash(st)
    fields = fieldnames(typeof(st))
    h = zero(UInt)
    for field in fields
        field_value = getfield(st, field)
        h = hash(field_value, h)
    end
    return h
end

function replace_nan!(matrix)
    rows, cols = size(matrix)
    distance = max(rows, cols)
    for i in 1:rows
        for j in 1:cols
            if isnan(matrix[i, j])
                neighbors = []
                for d in 1:distance
                    found = false
                    if i-d >= 1 && !isnan(matrix[i-d, j]);
                        push!(neighbors, matrix[i-d, j])
                        found = true
                    end
                    if i+d <= rows && !isnan(matrix[i+d, j])
                        push!(neighbors, matrix[i+d, j])
                        found = true
                    end
                    if j-d >= 1 && !isnan(matrix[i, j-d])
                        push!(neighbors, matrix[i, j-d])
                        found = true
                    end
                    if j+d <= cols && !isnan(matrix[i, j+d])
                        push!(neighbors, matrix[i, j+d])
                        found = true
                    end
                    if found; break; end
                end
                if !isempty(neighbors)
                    matrix[i, j] = sum(neighbors) / length(neighbors)
                end
            end
        end
    end
    return nothing
end